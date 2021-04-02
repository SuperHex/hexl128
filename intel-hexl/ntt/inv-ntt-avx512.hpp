// Copyright (C) 2020-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <immintrin.h>

#include <functional>
#include <vector>

#include "intel-hexl/ntt/ntt.hpp"
#include "logging/logging.hpp"
#include "ntt/ntt-avx512-util.hpp"
#include "ntt/ntt-internal.hpp"
#include "number-theory/number-theory.hpp"
#include "util/avx512-util.hpp"

namespace intel {
namespace hexl {

#ifdef HEXL_HAS_AVX512DQ

/// @brief The Harvey butterfly: assume X, Y in [0, 2p), and return X', Y' in
/// [0, 2p). such that X', Y' = X + Y (mod p), W(X - Y) (mod p).
/// @param[in,out] X Input representing 8 64-bit signed integers in SIMD form
/// @param[in,out] Y Input representing 8 64-bit signed integers in SIMD form
/// @param[in] W_op Root of unity representing 8 64-bit signed integers in SIMD
/// form
/// @param[in] W_precon Preconditioned \p W_op for BitShift-bit Barrett
/// reduction
/// @param[in] neg_modulus Negative modulus, i.e. (-p) represented as 8 64-bit
/// signed integers in SIMD form
/// @param[in] twice_modulus Twice the modulus, i.e. 2*p represented as 8 64-bit
/// signed integers in SIMD form
/// @param InputLessThanMod If true, assumes \p X, \p Y < \p p. Otherwise,
/// assumes \p X, \p Y < 2*\p p
/// @details See Algorithm 3 of https://arxiv.org/pdf/1205.2926.pdf
template <int BitShift, bool InputLessThanMod>
inline void InvButterfly(__m512i* X, __m512i* Y, __m512i W_op, __m512i W_precon,
                         __m512i neg_modulus, __m512i twice_modulus) {
  __m512i Y_minus_2p = _mm512_sub_epi64(*Y, twice_modulus);
  __m512i T = _mm512_sub_epi64(*X, Y_minus_2p);

  if (InputLessThanMod) {
    // No need for modulus reduction, since inputs are in [0,p)
    *X = _mm512_add_epi64(*X, *Y);
  } else {
    *X = _mm512_add_epi64(*X, Y_minus_2p);
    __mmask8 sign_bits = _mm512_movepi64_mask(*X);
    *X = _mm512_mask_add_epi64(*X, sign_bits, *X, twice_modulus);
  }
  __m512i Q = _mm512_hexl_mulhi_epi<BitShift>(W_precon, T);
  __m512i Q_p = _mm512_hexl_mullo_epi<BitShift>(Q, neg_modulus);
  *Y = _mm512_hexl_mullo_add_epi<BitShift>(Q_p, W_op, T);

  if (BitShift == 52) {
    // Discard high 12 bits; deals with case when W*T < Q*p in the low BitShift
    // bits.
    *Y = _mm512_and_epi64(*Y, _mm512_set1_epi64((1ULL << 52) - 1));
  }
}

template <int BitShift, bool InputLessThanMod>
void InvT1(uint64_t* operand, __m512i v_neg_modulus, __m512i v_twice_mod,
           uint64_t m, const uint64_t* W_op, const uint64_t* W_precon) {
  const __m512i* v_W_op_pt = reinterpret_cast<const __m512i*>(W_op);
  const __m512i* v_W_precon_pt = reinterpret_cast<const __m512i*>(W_precon);
  size_t j1 = 0;

  // 8 | m guaranteed by n >= 16
  HEXL_LOOP_UNROLL_8
  for (size_t i = m / 8; i > 0; --i) {
    uint64_t* X = operand + j1;
    __m512i* v_X_pt = reinterpret_cast<__m512i*>(X);

    __m512i v_X;
    __m512i v_Y;
    LoadInvInterleavedT1(X, &v_X, &v_Y);

    __m512i v_W_op = _mm512_loadu_si512(v_W_op_pt++);
    __m512i v_W_precon = _mm512_loadu_si512(v_W_precon_pt++);

    InvButterfly<BitShift, InputLessThanMod>(&v_X, &v_Y, v_W_op, v_W_precon,
                                             v_neg_modulus, v_twice_mod);

    _mm512_storeu_si512(v_X_pt++, v_X);
    _mm512_storeu_si512(v_X_pt, v_Y);

    j1 += 16;
  }
}

template <int BitShift>
void InvT2(uint64_t* X, __m512i v_neg_modulus, __m512i v_twice_mod, uint64_t m,
           const uint64_t* W_op, const uint64_t* W_precon) {
  // 4 | m guaranteed by n >= 16
  HEXL_LOOP_UNROLL_4
  for (size_t i = m / 4; i > 0; --i) {
    __m512i* v_X_pt = reinterpret_cast<__m512i*>(X);

    __m512i v_X;
    __m512i v_Y;
    LoadInvInterleavedT2(X, &v_X, &v_Y);

    __m512i v_W_op = LoadWOpT2(static_cast<const void*>(W_op));
    __m512i v_W_precon = LoadWOpT2(static_cast<const void*>(W_precon));

    InvButterfly<BitShift, false>(&v_X, &v_Y, v_W_op, v_W_precon, v_neg_modulus,
                                  v_twice_mod);

    _mm512_storeu_si512(v_X_pt++, v_X);
    _mm512_storeu_si512(v_X_pt, v_Y);
    X += 16;

    W_op += 4;
    W_precon += 4;
  }
}

template <int BitShift>
void InvT4(uint64_t* operand, __m512i v_neg_modulus, __m512i v_twice_mod,
           uint64_t m, const uint64_t* W_op, const uint64_t* W_precon) {
  uint64_t* X = operand;

  // 2 | m guaranteed by n >= 16
  HEXL_LOOP_UNROLL_4
  for (size_t i = m / 2; i > 0; --i) {
    __m512i* v_X_pt = reinterpret_cast<__m512i*>(X);

    __m512i v_X;
    __m512i v_Y;
    LoadInvInterleavedT4(X, &v_X, &v_Y);

    __m512i v_W_op = LoadWOpT4(static_cast<const void*>(W_op));
    __m512i v_W_precon = LoadWOpT4(static_cast<const void*>(W_precon));

    InvButterfly<BitShift, false>(&v_X, &v_Y, v_W_op, v_W_precon, v_neg_modulus,
                                  v_twice_mod);

    WriteInvInterleavedT4(v_X, v_Y, v_X_pt);
    X += 16;

    W_op += 2;
    W_precon += 2;
  }
}

template <int BitShift>
void InvT8(uint64_t* operand, __m512i v_neg_modulus, __m512i v_twice_mod,
           uint64_t t, uint64_t m, const uint64_t* W_op,
           const uint64_t* W_precon) {
  size_t j1 = 0;

  HEXL_LOOP_UNROLL_4
  for (size_t i = 0; i < m; i++) {
    uint64_t* X = operand + j1;
    uint64_t* Y = X + t;

    __m512i v_W_op = _mm512_set1_epi64(static_cast<int64_t>(*W_op++));
    __m512i v_W_precon = _mm512_set1_epi64(static_cast<int64_t>(*W_precon++));

    __m512i* v_X_pt = reinterpret_cast<__m512i*>(X);
    __m512i* v_Y_pt = reinterpret_cast<__m512i*>(Y);

    // assume 8 | t
    for (size_t j = t / 8; j > 0; --j) {
      __m512i v_X = _mm512_loadu_si512(v_X_pt);
      __m512i v_Y = _mm512_loadu_si512(v_Y_pt);

      InvButterfly<BitShift, false>(&v_X, &v_Y, v_W_op, v_W_precon,
                                    v_neg_modulus, v_twice_mod);

      _mm512_storeu_si512(v_X_pt++, v_X);
      _mm512_storeu_si512(v_Y_pt++, v_Y);
    }
    j1 += (t << 1);
  }
}

template <int BitShift>
void InverseTransformFromBitReverseAVX512(
    uint64_t* operand, uint64_t n, uint64_t mod,
    const uint64_t* inv_root_of_unity_powers,
    const uint64_t* precon_inv_root_of_unity_powers, uint64_t input_mod_factor,
    uint64_t output_mod_factor) {
  HEXL_CHECK(CheckNTTArguments(n, mod), "");
  HEXL_CHECK(mod < MaximumValue(BitShift) / 2,
             "mod " << mod << " too large for BitShift " << BitShift
                    << " => maximum value " << MaximumValue(BitShift) / 2);
  HEXL_CHECK_BOUNDS(precon_inv_root_of_unity_powers, n, MaximumValue(BitShift),
                    "precon_inv_root_of_unity_powers too large");
  HEXL_CHECK_BOUNDS(operand, n, MaximumValue(BitShift), "operand too large");
  HEXL_CHECK_BOUNDS(operand, n, input_mod_factor * mod,
                    "operand larger than input_mod_factor * modulus ("
                        << input_mod_factor << " * " << mod << ")");
  HEXL_CHECK(input_mod_factor == 1 || input_mod_factor == 2,
             "input_mod_factor must be 1 or 2; got " << input_mod_factor);
  HEXL_CHECK(output_mod_factor == 1 || output_mod_factor == 2,
             "output_mod_factor must be 1 or 2; got " << output_mod_factor);

  uint64_t twice_mod = mod << 1;
  __m512i v_modulus = _mm512_set1_epi64(static_cast<int64_t>(mod));
  __m512i v_neg_modulus = _mm512_set1_epi64(-static_cast<int64_t>(mod));
  __m512i v_twice_mod = _mm512_set1_epi64(static_cast<int64_t>(twice_mod));

  size_t t = 1;
  size_t root_index = 1;
  size_t m = (n >> 1);

  // Extract t=1, t=2, t=4 loops separately
  {
    // t = 1
    const uint64_t* W_op = &inv_root_of_unity_powers[root_index];
    const uint64_t* W_precon = &precon_inv_root_of_unity_powers[root_index];
    if (input_mod_factor == 1) {
      InvT1<BitShift, true>(operand, v_neg_modulus, v_twice_mod, m, W_op,
                            W_precon);
    } else {
      InvT1<BitShift, false>(operand, v_neg_modulus, v_twice_mod, m, W_op,
                             W_precon);
    }
    t <<= 1;
    root_index += m;
    m >>= 1;

    // t = 2
    W_op = &inv_root_of_unity_powers[root_index];
    W_precon = &precon_inv_root_of_unity_powers[root_index];
    InvT2<BitShift>(operand, v_neg_modulus, v_twice_mod, m, W_op, W_precon);

    t <<= 1;
    root_index += m;
    m >>= 1;

    // t = 4
    W_op = &inv_root_of_unity_powers[root_index];
    W_precon = &precon_inv_root_of_unity_powers[root_index];
    InvT4<BitShift>(operand, v_neg_modulus, v_twice_mod, m, W_op, W_precon);
    t <<= 1;
    root_index += m;
    m >>= 1;
  }

  // t >= 8
  for (; m > 1; m >>= 1) {
    const uint64_t* W_op = &inv_root_of_unity_powers[root_index];
    const uint64_t* W_precon = &precon_inv_root_of_unity_powers[root_index];
    InvT8<BitShift>(operand, v_neg_modulus, v_twice_mod, t, m, W_op, W_precon);
    t <<= 1;
    root_index += m;
  }

  HEXL_VLOG(4, "AVX512 intermediate operand "
                   << std::vector<uint64_t>(operand, operand + n));

  const uint64_t W_op = inv_root_of_unity_powers[root_index];
  MultiplyFactor mf_inv_n(InverseUIntMod(n, mod), BitShift, mod);
  const uint64_t inv_n = mf_inv_n.Operand();
  const uint64_t inv_n_prime = mf_inv_n.BarrettFactor();

  MultiplyFactor mf_inv_n_w(MultiplyUIntMod(inv_n, W_op, mod), BitShift, mod);
  const uint64_t inv_n_w = mf_inv_n_w.Operand();
  const uint64_t inv_n_w_prime = mf_inv_n_w.BarrettFactor();

  HEXL_VLOG(4, "inv_n_w " << inv_n_w);

  uint64_t* X = operand;
  uint64_t* Y = X + (n >> 1);

  __m512i v_inv_n = _mm512_set1_epi64(static_cast<int64_t>(inv_n));
  __m512i v_inv_n_prime = _mm512_set1_epi64(static_cast<int64_t>(inv_n_prime));
  __m512i v_inv_n_w = _mm512_set1_epi64(static_cast<int64_t>(inv_n_w));
  __m512i v_inv_n_w_prime =
      _mm512_set1_epi64(static_cast<int64_t>(inv_n_w_prime));

  __m512i* v_X_pt = reinterpret_cast<__m512i*>(X);
  __m512i* v_Y_pt = reinterpret_cast<__m512i*>(Y);

  const __m512i two_pow52_min1 = _mm512_set1_epi64((1ULL << 52) - 1);

  // Merge final InvNTT loop with modulus reduction baked-in
  HEXL_LOOP_UNROLL_4
  for (size_t j = n / 16; j > 0; --j) {
    __m512i v_X = _mm512_loadu_si512(v_X_pt);
    __m512i v_Y = _mm512_loadu_si512(v_Y_pt);

    // Slightly different from regular InvButterfly because different W is used
    // for X and Y

    __m512i Y_minus_2p = _mm512_sub_epi64(v_Y, v_twice_mod);
    __m512i X_plus_Y_mod2p =
        _mm512_hexl_small_add_mod_epi64(v_X, v_Y, v_twice_mod);
    // T = *X + twice_mod - *Y
    __m512i T = _mm512_sub_epi64(v_X, Y_minus_2p);

    __m512i Q1 = _mm512_hexl_mulhi_epi<BitShift>(v_inv_n_prime, X_plus_Y_mod2p);
    // X = inv_N * X_plus_Y_mod2p - Q1 * modulus;
    __m512i inv_N_tx = _mm512_hexl_mullo_epi<BitShift>(v_inv_n, X_plus_Y_mod2p);
    v_X = _mm512_hexl_mullo_add_epi<BitShift>(inv_N_tx, Q1, v_neg_modulus);
    if (BitShift == 52) {
      // Discard high 12 bits; deals with case when W*T < Q1*p in the low
      // BitShift bits.
      v_X = _mm512_and_epi64(v_X, two_pow52_min1);
    }

    __m512i Q2 = _mm512_hexl_mulhi_epi<BitShift>(v_inv_n_w_prime, T);
    // Y = inv_N_W * T - Q2 * modulus;
    __m512i inv_N_W_T = _mm512_hexl_mullo_epi<BitShift>(v_inv_n_w, T);
    v_Y = _mm512_hexl_mullo_add_epi<BitShift>(inv_N_W_T, Q2, v_neg_modulus);
    if (BitShift == 52) {
      // Discard high 12 bits; deals with case when W*T < Q2*p in the low
      // BitShift bits.
      v_Y = _mm512_and_epi64(v_Y, two_pow52_min1);
    }

    if (output_mod_factor == 1) {
      // Modulus reduction from [0,2p), to [0,p)
      v_X = _mm512_hexl_small_mod_epu64(v_X, v_modulus);
      v_Y = _mm512_hexl_small_mod_epu64(v_Y, v_modulus);
    }

    _mm512_storeu_si512(v_X_pt++, v_X);
    _mm512_storeu_si512(v_Y_pt++, v_Y);
  }

  HEXL_VLOG(5, "AVX512 returning operand "
                   << std::vector<uint64_t>(operand, operand + n));
}

#endif  // HEXL_HAS_AVX512DQ

}  // namespace hexl
}  // namespace intel

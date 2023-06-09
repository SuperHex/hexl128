// Copyright (C) 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cmath>

#include "hexl/util/check.hpp"
#include "hexl/util/types.hpp"

namespace intel {
namespace hexl {

#ifdef HEXL_USE_GNU
// Return x * y as 128-bit integer
// Correctness if x * y < 128 bits
inline uint256_t MultiplyUInt128(uint128_t x, uint128_t y) {
  return uint256_t(x) * uint256_t(y);
}

inline uint128_t BarrettReduce256(uint128_t input_hi, uint128_t input_lo,
                                 uint128_t modulus) {
  HEXL_CHECK(modulus != 0, "modulus == 0")
  uint256_t n = (static_cast<uint256_t>(input_hi) << 128) |
                (static_cast<uint256_t>(input_lo));

  return static_cast<uint128_t>(n % modulus);
  // TODO(fboemer): actually use barrett reduction if performance-critical
}

// Returns low 64bit of 128b/64b where x1=high 64b, x0=low 64b
inline uint128_t DivideUInt256UInt128Lo(uint128_t x1, uint128_t x0, uint128_t y) {
  uint256_t n =
      (static_cast<uint256_t>(x1) << 128) | (static_cast<uint256_t>(x0));
  uint256_t q = n / y;

  return static_cast<uint128_t>(q);
}

// Multiplies x * y as 128-bit integer.
// @param prod_hi Stores high 64 bits of product
// @param prod_lo Stores low 64 bits of product
inline void MultiplyUInt128(uint128_t x, uint128_t y, uint128_t* prod_hi,
                           uint128_t* prod_lo) {
  uint256_t prod = MultiplyUInt128(x, y);
  *prod_hi = static_cast<uint128_t>(prod >> 128);
  *prod_lo = static_cast<uint128_t>(prod);
}

// Return the high 128 minus BitShift bits of the 128-bit product x * y
template <int BitShift>
inline uint128_t MultiplyUInt128Hi(uint128_t x, uint128_t y) {
  uint256_t product = MultiplyUInt128(x, y);
  return static_cast<uint128_t>(product >> BitShift);
}

// Returns most-significant bit of the input
inline uint128_t MSB(uint128_t input) {
    return static_cast<uint128_t>(mp::msb(input));
}

#define HEXL_LOOP_UNROLL_4 _Pragma("GCC unroll 4")
#define HEXL_LOOP_UNROLL_8 _Pragma("GCC unroll 8")

#endif

}  // namespace hexl
}  // namespace intel

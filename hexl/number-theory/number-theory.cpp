// Copyright (C) 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#include "hexl/number-theory/number-theory.hpp"

#include "hexl/logging/logging.hpp"
#include "hexl/util/check.hpp"
#include "util/util-internal.hpp"

namespace intel {
namespace hexl {

uint128_t InverseMod(uint128_t input, uint128_t modulus) {
  uint128_t a = input % modulus;
  HEXL_CHECK(a != 0, input << " does not have a InverseMod");

  if (modulus == 1) {
    return 0;
  }

  int128_t m0 = static_cast<int128_t>(modulus);
  int128_t y = 0;
  int128_t x = 1;
  while (a > 1) {
    // q is quotient
    int128_t q = static_cast<int128_t>(a / modulus);

    int128_t t = static_cast<int128_t>(modulus);
    modulus = a % modulus;
    a = static_cast<uint128_t>(t);

    // Update y and x
    t = y;
    y = x - q * y;
    x = t;
  }

  // Make x positive
  if (x < 0) x += m0;

  return uint128_t(x);
}

uint128_t MultiplyMod(uint128_t x, uint128_t y, uint128_t modulus) {
  HEXL_CHECK(modulus != 0, "modulus == 0");
  HEXL_CHECK(x < modulus, "x " << x << " >= modulus " << modulus);
  HEXL_CHECK(y < modulus, "y " << y << " >= modulus " << modulus);
  uint128_t prod_hi, prod_lo;
  MultiplyUInt128(x, y, &prod_hi, &prod_lo);

  return BarrettReduce256(prod_hi, prod_lo, modulus);
}

uint128_t MultiplyMod(uint128_t x, uint128_t y, uint128_t y_precon,
                     uint128_t modulus) {
  uint128_t q = MultiplyUInt128Hi<128>(x, y_precon);
  q = x * y - q * modulus;
  return q >= modulus ? q - modulus : q;
}

uint128_t AddUIntMod(uint128_t x, uint128_t y, uint128_t modulus) {
  HEXL_CHECK(x < modulus, "x " << x << " >= modulus " << modulus);
  HEXL_CHECK(y < modulus, "y " << y << " >= modulus " << modulus);
  uint128_t sum = x + y;
  return (sum >= modulus) ? (sum - modulus) : sum;
}

uint128_t SubUIntMod(uint128_t x, uint128_t y, uint128_t modulus) {
  HEXL_CHECK(x < modulus, "x " << x << " >= modulus " << modulus);
  HEXL_CHECK(y < modulus, "y " << y << " >= modulus " << modulus);
  uint128_t diff = (x + modulus) - y;
  return (diff >= modulus) ? (diff - modulus) : diff;
}

// Returns base^exp mod modulus
uint128_t PowMod(uint128_t base, uint128_t exp, uint128_t modulus) {
  base %= modulus;
  uint128_t result = 1;
  while (exp > 0) {
    if (exp & 1) {
      result = MultiplyMod(result, base, modulus);
    }
    base = MultiplyMod(base, base, modulus);
    exp >>= 1;
  }
  return result;
}

// Returns true whether root is a degree-th root of unity
// degree must be a power of two.
bool IsPrimitiveRoot(uint128_t root, uint128_t degree, uint128_t modulus) {
  if (root == 0) {
    return false;
  }
  HEXL_CHECK(IsPowerOfTwo(degree), degree << " not a power of 2");

  HEXL_VLOG(4, "IsPrimitiveRoot root " << root << ", degree " << degree
                                       << ", modulus " << modulus);

  // Check if root^(degree/2) == -1 mod modulus
  return PowMod(root, degree / 2, modulus) == (modulus - 1);
}

// Tries to return a primitive degree-th root of unity
// throw error if no root is found
uint128_t GeneratePrimitiveRoot(uint128_t degree, uint128_t modulus) {
  // We need to divide modulus-1 by degree to get the size of the quotient group
  uint128_t size_entire_group = modulus - 1;

  // Compute size of quotient group
  uint128_t size_quotient_group = size_entire_group / degree;

  for (int trial = 0; trial < 200; ++trial) {
    uint128_t root = GenerateInsecureUniformIntRandomValue(0, modulus);
    root = PowMod(root, size_quotient_group, modulus);

    if (IsPrimitiveRoot(root, degree, modulus)) {
      return root;
    }
  }
  HEXL_CHECK(false, "no primitive root found for degree "
                        << degree << " modulus " << modulus);
  return 0;
}

// Returns true whether root is a degree-th root of unity
// degree must be a power of two.
uint128_t MinimalPrimitiveRoot(uint128_t degree, uint128_t modulus) {
  HEXL_CHECK(IsPowerOfTwo(degree),
             "Degere " << degree << " is not a power of 2");

  uint128_t root = GeneratePrimitiveRoot(degree, modulus);

  uint128_t generator_sq = MultiplyMod(root, root, modulus);
  uint128_t current_generator = root;

  uint128_t min_root = root;

  // Check if root^(degree/2) == -1 mod modulus
  for (size_t i = 0; i < degree; ++i) {
    if (current_generator < min_root) {
      min_root = current_generator;
    }
    current_generator = MultiplyMod(current_generator, generator_sq, modulus);
  }

  return min_root;
}

uint64_t ReverseBits(uint64_t x, uint64_t bit_width) {
  HEXL_CHECK(x == 0 || MSB(x) <= bit_width, "MSB(" << x << ") = " << MSB(x)
                                                   << " must be >= bit_width "
                                                   << bit_width)
  if (bit_width == 0) {
    return 0;
  }
  uint64_t rev = 0;
  for (uint64_t i = bit_width; i > 0; i--) {
    rev |= ((x & 1) << (i - 1));
    x >>= 1;
  }
  return rev;
}

// Miller-Rabin primality test
bool IsPrime(uint128_t n) {
  // n < 2^64, so it is enough to test a=2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
  // and 37. See
  // https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
  static const std::vector<uint128_t> as{2,  3,  5,  7,  11, 13,
                                        17, 19, 23, 29, 31, 37};

  for (const uint128_t& a : as) {
    if (n == a) return true;
    if (n % a == 0) return false;
  }

  // Write n == 2**r * d + 1 with d odd.
  size_t r = 127;
  while (r > 0) {
      uint128_t two_pow_r = (uint128_t{1ULL} << r);
    if ((n - 1) % two_pow_r == 0) {
      break;
    }
    --r;
  }
  HEXL_CHECK(r != 0, "Error factoring n " << n);
  uint128_t d = (n - 1) / (uint128_t{1ULL} << r);

  HEXL_CHECK(n == (uint128_t{1ULL} << r) * d + 1, "Error factoring n " << n);
  HEXL_CHECK(d % 2 == 1, "d is even");

  for (const uint128_t& a : as) {
    uint128_t x = PowMod(a, d, n);
    if ((x == 1) || (x == n - 1)) {
      continue;
    }

    bool prime = false;
    for (uint128_t i = 1; i < r; ++i) {
      x = PowMod(x, 2, n);
      if (x == n - 1) {
        prime = true;
        break;
      }
    }
    if (!prime) {
      return false;
    }
  }
  return true;
}

std::vector<uint128_t> GeneratePrimes(size_t num_primes, size_t bit_size,
                                     bool prefer_small_primes,
                                     size_t ntt_size) {
  HEXL_CHECK(num_primes > 0, "num_primes == 0");
  HEXL_CHECK(IsPowerOfTwo(ntt_size),
             "ntt_size " << ntt_size << " is not a power of two");
  HEXL_CHECK(Log2(ntt_size) < bit_size,
             "log2(ntt_size) " << Log2(ntt_size)
                               << " should be less than bit_size " << bit_size);

  int128_t prime_lower_bound = (int128_t{1LL} << bit_size) + 1LL;
  int128_t prime_upper_bound = (int128_t{1LL} << (bit_size + 1LL)) - 1LL;

  // Keep signed to enable negative step
  int128_t prime_candidate =
      prefer_small_primes
          ? prime_lower_bound
          : prime_upper_bound - (prime_upper_bound % (2 * ntt_size)) + 1;
  HEXL_CHECK(prime_candidate % (2 * ntt_size) == 1, "bad prime candidate");

  // Ensure prime % 2 * ntt_size == 1
  int128_t prime_candidate_step =
      (prefer_small_primes ? 1 : -1) * 2 * static_cast<int128_t>(ntt_size);

  auto continue_condition = [&](int128_t local_candidate_prime) {
    if (prefer_small_primes) {
      return local_candidate_prime < prime_upper_bound;
    } else {
      return local_candidate_prime > prime_lower_bound;
    }
  };

  std::vector<uint128_t> ret;

  while (continue_condition(prime_candidate)) {
      if (IsPrime(static_cast<uint128_t>(prime_candidate))) {
      HEXL_CHECK(prime_candidate % (2 * ntt_size) == 1, "bad prime candidate");
      ret.emplace_back(static_cast<uint128_t>(prime_candidate));
      if (ret.size() == num_primes) {
        return ret;
      }
    }
    prime_candidate += prime_candidate_step;
  }

  HEXL_CHECK(false, "Failed to find enough primes");
  return ret;
}

}  // namespace hexl
}  // namespace intel

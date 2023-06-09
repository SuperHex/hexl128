// Copyright (C) 2020 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdint.h>

#include "hexl/util/defines.hpp"

// #if defined(HEXL_USE_GNU) || defined(HEXL_USE_CLANG)
// __extension__ typedef __int128 int128_t;
// __extension__ typedef unsigned __int128 uint128_t;
// #endif

#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

using int128_t  = mp::int128_t;
using uint128_t = mp::uint128_t;
using int256_t  = mp::int256_t;
using uint256_t = mp::uint256_t;


/*
 * Copyright 2020 Improbable Worlds Limited
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef PHTREE_COMMON_COMMON_H
#define PHTREE_COMMON_COMMON_H

#include "b_plus_tree_map.h"
#include "base_types.h"
#include "bits.h"
#include "flat_array_map.h"
#include "flat_sparse_map.h"
#include "tree_stats.h"
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>

namespace improbable::phtree {

// This is the single-point inclusion file for common types/function/... for the PH-Tree.
// 'single-point inclusion' meaning that including it provides all relevant types/functions/... .

// ************************************************************************
// Bits
// ************************************************************************

/*
 * Encode the bits at the given position of all attributes into a hyper-cube address.
 * Currently, the first attribute determines the left-most (high-value) bit of the address
 * (left to right ordered)
 *
 * @param valSet vector
 * @param postfix_len the postfix length
 * @returns Encoded HC position, which is the index in the array if the entries would be stored in
 * an array.
 */
template <dimension_t DIM, typename SCALAR>
static hc_pos_t CalcPosInArray(const PhPoint<DIM, SCALAR>& valSet, bit_width_t postfix_len) {
    // n=DIM,  i={0..n-1}
    // i = 0 :  |0|1|0|1|0|1|0|1|
    // i = 1 :  | 0 | 1 | 0 | 1 |
    // i = 2 :  |   0   |   1   |
    // len = 2^n
    // Following formula was for inverse ordering of current ordering...
    // pos = sum (i=1..n, len/2^i) = sum (..., 2^(n-i))
    bit_mask_t<SCALAR> valMask = bit_mask_t<SCALAR>(1) << postfix_len;
    hc_pos_t pos = 0;
    for (dimension_t i = 0; i < DIM; ++i) {
        pos <<= 1;
        // set pos-bit if bit is set in value
        pos |= (valMask & valSet[i]) >> postfix_len;
    }
    return pos;
}

template <dimension_t DIM, typename SCALAR>
static bool IsInRange(
    const PhPoint<DIM, SCALAR>& candidate,
    const PhPoint<DIM, SCALAR>& range_min,
    const PhPoint<DIM, SCALAR>& range_max) {
    for (dimension_t i = 0; i < DIM; ++i) {
        auto k = candidate[i];
        if (k < range_min[i] || k > range_max[i]) {
            return false;
        }
    }
    return true;
}

/*
 * @param v1 key 1
 * @param v2 key 2
 * @return the number of diverging bits. For each dimension we determine the most significant bit
 * where the two keys differ. We then count this bit plus all trailing bits (even if individual bits
 * may be the same). Then we return the highest number of diverging bits found in any dimension of
 * the two keys. In case of key1==key2 we return 0. In other words, for 64 bit keys, we return 64
 * minus the number of leading bits that are common in both keys across all dimensions.
 */
template <dimension_t DIM, typename SCALAR>
static bit_width_t NumberOfDivergingBits(
    const PhPoint<DIM, SCALAR>& v1, const PhPoint<DIM, SCALAR>& v2) {
    // write all differences to diff, we just check diff afterwards
    bit_mask_t<SCALAR> diff = 0;
    for (dimension_t i = 0; i < DIM; ++i) {
        diff |= (v1[i] ^ v2[i]);
    }
    assert(CountLeadingZeros(diff) <= MAX_BIT_WIDTH<SCALAR>);
    return MAX_BIT_WIDTH<SCALAR> - CountLeadingZeros(diff);
}

template <dimension_t DIM, typename SCALAR>
static bool KeyEquals(
    const PhPoint<DIM, SCALAR>& key_a, const PhPoint<DIM, SCALAR>& key_b, bit_mask_t<SCALAR> mask) {
    for (dimension_t i = 0; i < DIM; ++i) {
        if (((key_a[i] ^ key_b[i]) & mask) != 0) {
            return false;
        }
    }
    return true;
}

// ************************************************************************
// String helpers
// ************************************************************************

template <typename SCALAR>
static inline std::string ToBinary(SCALAR l, bit_width_t width = MAX_BIT_WIDTH<SCALAR>) {
    std::ostringstream sb;
    // long mask = DEPTH < 64 ? (1<<(DEPTH-1)) : 0x8000000000000000L;
    for (bit_width_t i = 0; i < width; ++i) {
        bit_mask_t<SCALAR> mask = (bit_mask_t<SCALAR>(1) << (width - i - 1));
        sb << ((l & mask) != 0 ? "1" : "0");
        if ((i + 1) % 8 == 0 && (i + 1) < width) {
            sb << '.';
        }
    }
    return sb.str();
}

template <dimension_t DIM, typename SCALAR>
static inline std::string ToBinary(
    const PhPoint<DIM, SCALAR>& la, bit_width_t width = MAX_BIT_WIDTH<SCALAR>) {
    std::ostringstream sb;
    for (dimension_t i = 0; i < DIM; ++i) {
        sb << ToBinary(la[i], width) << ", ";
    }
    return sb.str();
}

}  // namespace improbable::phtree

#endif  // PHTREE_COMMON_COMMON_H

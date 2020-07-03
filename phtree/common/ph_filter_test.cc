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

#include "ph_common.h"
#include <gtest/gtest.h>
#include <random>

using namespace improbable::phtree;

TEST(PhTreeFilterTest, BoxFilterTest) {
    auto filter = PhFilterAABB<2, int*>({3, 3}, {7, 7});
    // root is always valid
    ASSERT_TRUE(filter.IsNodeValid({0, 0}, 63));
    // valid because node encompasses the AABB
    ASSERT_TRUE(filter.IsNodeValid({1, 1}, 10));
    // valid
    ASSERT_TRUE(filter.IsNodeValid({7, 7}, 1));
    // invalid
    ASSERT_FALSE(filter.IsNodeValid({88, 5}, 1));

    ASSERT_TRUE(filter.IsEntryValid({3, 7}, nullptr));
    ASSERT_FALSE(filter.IsEntryValid({2, 8}, nullptr));
}

TEST(PhTreeFilterTest, FilterNoOpSmokeTest) {
    auto filter = PhFilterNoOp<3, int>();
    ASSERT_TRUE(filter.IsNodeValid({3, 7, 2}, 10));
    ASSERT_TRUE(filter.IsEntryValid({3, 7, 2}, 10));
}
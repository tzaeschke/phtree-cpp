/*
 * Copyright 2022 Tilmann Zäschke
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

#include "b_plus_tree_hash_map.h"
#include <gtest/gtest.h>
#include <random>
#include <unordered_set>

using namespace improbable::phtree;

static int default_construct_count_ = 0;
static int construct_count_ = 0;
static int copy_construct_count_ = 0;
static int move_construct_count_ = 0;
static int copy_assign_count_ = 0;
static int move_assign_count_ = 0;
static int destruct_count_ = 0;

static void reset_id_counters() {
    default_construct_count_ = 0;
    construct_count_ = 0;
    copy_construct_count_ = 0;
    move_construct_count_ = 0;
    copy_assign_count_ = 0;
    move_assign_count_ = 0;
    destruct_count_ = 0;
}

static void print_id_counters() {
    std::cout << "dc=" << default_construct_count_ << " c=" << construct_count_
              << " cc=" << copy_construct_count_ << " mc=" << move_construct_count_
              << " ca=" << copy_assign_count_ << " ma=" << move_assign_count_
              << " d=" << destruct_count_ << std::endl;
}

struct Id {
    Id() : _i{0} {
        ++default_construct_count_;
    }

    explicit Id(const size_t i) : _i{static_cast<int>(i)} {
        ++construct_count_;
    }

    explicit Id(const int i) : _i{i} {
        ++construct_count_;
    }

    Id(const Id& other) {
        ++copy_construct_count_;
        _i = other._i;
    }

    Id(Id&& other) noexcept {
        ++move_construct_count_;
        _i = other._i;
    }

    Id& operator=(const Id& other) noexcept {
        ++copy_assign_count_;
        _i = other._i;
        return *this;
    }
    Id& operator=(Id&& other) noexcept {
        ++move_assign_count_;
        _i = other._i;
        return *this;
    }

    bool operator==(const Id& rhs) const {
        return _i == rhs._i;
    }

    bool operator==(Id&& rhs) const {
        return _i == rhs._i;
    }

    ~Id() {
        ++destruct_count_;
    }

    int _i;
};

namespace std {
template <>
struct hash<Id> {
    size_t operator()(const Id& x) const {
        return std::hash<int>{}(x._i % 10);
    }
};
};  // namespace std

//struct IdHash {
//    template <class T1, class T2>
//    std::size_t operator()(std::pair<T1, T2> const& v) const {
//        return std::hash<T1>()(v.size());
//    }
//};

TEST(PhTreeBptHashMapTest, SmokeTest) {
    const int max_size = 200;

    std::default_random_engine random_engine{0};
    std::uniform_int_distribution<> cube_distribution(0, max_size - 1);

    int val = 0;
    for (int i = 0; i < 10; i++) {
        b_plus_tree_hash_map<size_t, size_t, std::hash<size_t>, std::equal_to<size_t>, max_size> test_map;
        std::map<size_t, size_t> reference_map;
        for (int j = 0; j < 2 * max_size; j++) {
            size_t key = cube_distribution(random_engine);
            bool hasVal = test_map.find(key) != test_map.end();
            bool hasValRef = reference_map.find(key) != reference_map.end();
            ASSERT_EQ(hasVal, hasValRef);

            reference_map.emplace(key, val);
            test_map.emplace(key, val);
            test_map._check();

            std::cout << "i=" << i << " j=" << j << " k/v=" << key << "/" << val << std::endl;

            ASSERT_EQ(test_map.size(), reference_map.size());
            for (auto it : reference_map) {
                size_t kRef = it.first;
                size_t vMap = test_map.find(kRef)->second;
                ASSERT_EQ(vMap, it.second);
            }
            for (auto it : test_map) {
                size_t k = it.first;
                size_t vRef = reference_map.find(k)->second;
                size_t vMap = test_map.find(k)->second;
                ASSERT_EQ(vMap, vRef);
            }
            ++val;
        }
    }
}

TEST(PhTreeBptHashSetTest, SmokeTestNonUnique) {
    const int max_size = 200;

    std::default_random_engine random_engine{0};
    std::uniform_int_distribution<> cube_distribution(0, max_size - 1);

    for (int i = 0; i < 10; i++) {
        b_plus_tree_hash_set<Id> test_map;
        std::unordered_set<Id> reference_map;
        int n = 0;
        for (int j = 0; j < 2 * max_size; j++) {
            size_t key = cube_distribution(random_engine);
            Id id(key);
            std::cout << "i=" << i << " j=" << j << " k/n=" << key << "/" << n << std::endl;
            bool hasVal = test_map.find(id) != test_map.end();
            bool hasValRef = reference_map.find(id) != reference_map.end();
            ASSERT_EQ(hasVal, hasValRef);

            reference_map.emplace(id);
            test_map.emplace(id);
            test_map._check();

            ASSERT_EQ(test_map.size(), reference_map.size());
            for (auto id : reference_map) {
                Id& idMap = *test_map.find(id);
                ASSERT_EQ(idMap, id);
            }
            for (auto id : test_map) {
                const Id& vRef = *reference_map.find(id);
                Id& vMap = *test_map.find(id);  // TODO require const?
                ASSERT_EQ(vMap, vRef);
            }
            ++n;
        }
    }
}

TEST(PhTreeBptHashMapTest, SmokeTestWithTryEmplace) {
    const int max_size = 200;

    std::default_random_engine random_engine{0};
    std::uniform_int_distribution<> cube_distribution(0, max_size - 1);

    for (int i = 0; i < 10; i++) {
        b_plus_tree_hash_map<size_t, size_t, std::hash<size_t>, std::equal_to<size_t>, max_size> test_map;
        std::map<size_t, size_t> reference_map;
        for (int j = 0; j < 2 * max_size; j++) {
            size_t val = cube_distribution(random_engine);
            bool hasVal = test_map.find(val) != test_map.end();
            bool hasValRef = reference_map.find(val) != reference_map.end();
            ASSERT_EQ(hasVal, hasValRef);
            if (!hasVal) {
                reference_map.emplace(val, val);
                test_map.try_emplace(val, val);
            }
            ASSERT_EQ(test_map.size(), reference_map.size());
            for (auto it : reference_map) {
                size_t vRef = it.first;
                size_t vMap = test_map.find(vRef)->second;
                ASSERT_EQ(vMap, vRef);
            }
            for (auto it : test_map) {
                size_t v = it.first;
                size_t vRef = reference_map.find(v)->second;
                size_t vMap = test_map.find(v)->second;
                ASSERT_EQ(vMap, vRef);
            }
        }
    }
}

TEST(PhTreeBptHashMapTest, SmokeTestWithErase) {
    const int max_size = 200;

    std::default_random_engine random_engine{0};
    std::uniform_int_distribution<> cube_distribution(0, max_size - 1);

    for (int i = 0; i < 10; i++) {
        b_plus_tree_hash_map<size_t, Id, std::hash<size_t>, std::equal_to<size_t>, max_size> test_map{};
        std::unordered_map<size_t, size_t> reference_map{};
        std::vector<size_t> key_list{};
        for (int j = 0; j < 2 * max_size; j++) {
            size_t val = cube_distribution(random_engine);
            bool hasVal = test_map.find(val) != test_map.end();
            bool hasValRef = reference_map.find(val) != reference_map.end();
            ASSERT_EQ(hasVal, hasValRef);
            if (!hasVal) {
                reference_map.emplace(val, val);
                test_map.try_emplace(val, val);
                key_list.emplace_back(val);
            }
        }

        std::shuffle(key_list.begin(), key_list.end(), random_engine);
        for (auto key : key_list) {
            if (key % 2 == 0) {
                test_map.erase(key);
            } else {
                auto it = test_map.find(key);
                ASSERT_NE(it, test_map.end());
                ASSERT_EQ(it->second, key);
                test_map.erase(it);
            }
            test_map._check();
            reference_map.erase(key);
            for (auto it : reference_map) {
                size_t vRef = it.first;
                size_t vMap = test_map.find(vRef)->second;
                ASSERT_EQ(vMap, vRef);
            }
            for (auto it : test_map) {
                size_t v = it.first;
                size_t vRef = reference_map.find(v)->second;
                size_t vMap = test_map.find(v)->second;
                ASSERT_EQ(vMap, vRef);
            }
            ASSERT_EQ(test_map.size(), reference_map.size());
        }
    }
}

TEST(PhTreeBptHashMapTest, SmokeTestLowerBound) {
    const int max_size = 200;

    std::default_random_engine random_engine{0};
    std::uniform_int_distribution<> cube_distribution(0, max_size - 1);

    for (int i = 0; i < 10; i++) {
        b_plus_tree_hash_map<size_t, Id, std::hash<size_t>, std::equal_to<size_t>, max_size> test_map;
        std::map<size_t, size_t> reference_map;
        for (int j = 0; j < 2 * max_size; j++) {
            size_t val = cube_distribution(random_engine);
            bool hasVal = test_map.find(val) != test_map.end();
            bool hasValRef = reference_map.find(val) != reference_map.end();
            ASSERT_EQ(hasVal, hasValRef);
            if (!hasVal) {
                reference_map.emplace(val, val);
                test_map.try_emplace(val, val);
            }
            ASSERT_EQ(test_map.size(), reference_map.size());
            for (auto it : reference_map) {
                size_t vRef = it.first;
                size_t vMap = test_map.lower_bound(vRef)->second;
                ASSERT_EQ(vMap, vRef);
            }
            for (auto it : test_map) {
                size_t v = it.first;
                size_t vRef = reference_map.find(v)->second;
                size_t vMap = test_map.lower_bound(v)->second;
                ASSERT_EQ(vMap, vRef);
            }
            for (size_t v = 0; v < max_size + 5; ++v) {
                auto itRef = reference_map.lower_bound(v);
                auto itMap = test_map.lower_bound(v);
                if (itRef == reference_map.end()) {
                    ASSERT_EQ(itMap, test_map.end());
                } else {
                    ASSERT_NE(itMap, test_map.end());
                    // ASSERT_EQ(v, itRef->second);
                    ASSERT_EQ(itRef->second, itMap->second);
                }
            }
        }
    }
}

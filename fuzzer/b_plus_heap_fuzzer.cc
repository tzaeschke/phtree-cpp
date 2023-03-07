/*
 * Copyright 2023 Tilmann Zäschke
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

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <ostream>

#include <map>

#include "include/phtree/common/b_plus_tree_heap.h"

static volatile int Sink;

using Instruction = std::uint8_t;
using Key = std::uint8_t;
using Value = std::uint8_t;

constexpr bool PRINT = !true;

void print() {}

using Id = std::pair<uint8_t, uint8_t>;

struct CompareId {
    bool operator()(const Id& left, const Id& right) const {
        return left.first > right.first;
    };
};


extern "C" int LLVMFuzzerTestOneInput(const uint8_t* Data, size_t Size) {
    assert(Data);

    if (PRINT) {
        std::cout << "TEST(PhTreeBptHeapTest, FuzzTest1) {" << std::endl;
        std::cout << "    using Key = std::uint8_t;" << std::endl;
        std::cout << "    using Value = std::uint8_t;" << std::endl;
        std::cout << "    using Id = std::pair<uint8_t, uint8_t>;" << std::endl;
        std::cout << "    b_plus_tree_heap<Id, std::greater<double>> tree{};" << std::endl;
    }

    auto scopeguard = []() { std::cout << "};" << std::endl; };

    //phtree::bptree::b_plus_tree_heap<Id, CompareId> tree;
    phtree::bptree::b_plus_tree_heap<Id, std::greater<double>> tree;
    std::multimap<Key, Value, std::greater<double>> map;

    size_t pos = 0;

    while (pos + 4 < Size) {
        Instruction inst = Data[pos++] % 5;
        Key key = Data[pos++];
        Value value = Data[pos++];
        Id entry{key, value};
        switch (inst) {
        case 0: {
            if (PRINT)
                std::cout << "    tree.emplace(" << (int)key << ", " << (int)value << ");"
                          << std::endl;
            tree.emplace(entry);
            map.emplace(key, value);
            break;
        }
        case 1: {
            if (!tree.empty()) {
                if (PRINT)
                    std::cout << "    tree.pop();" << std::endl;
                tree.pop();
                map.erase(--map.end());
            }
            break;
        }
        case 2: {
            if (!tree.empty()) {
                if (PRINT)
                    std::cout << "    tree.pop_max();" << std::endl;
                tree.pop_max();
                map.erase(map.begin());
            }
            break;
        }
        case 3: {
            if (!tree.empty()) {
                if (PRINT)
                    std::cout << "    auto x = tree.top();" << std::endl;
                auto& x = tree.top();
                auto& x2 = *(--map.end());
                assert(x.first == x2.first);
            }
            break;
        }
        case 4: {
            if (!tree.empty()) {
                if (PRINT)
                    std::cout << "    auto x = tree.top_max();" << std::endl;
                auto& x = tree.top_max();
                auto& x2 = *(map.begin());
                assert(x.first == x2.first);
            }
            break;
        }
        default:
            std::cout << "Unexpected instruction: " << inst << std::endl;
        }
    }

    tree._check();

//    for (auto& entry : map) {
//        const Key& vRef = entry.first;
//        Key vMap = tree.find(vRef)->first;
//        assert(vMap == vRef);
//    }
//    for (auto& entry : tree) {
//        Key v = entry.first;
//        const Key& vRef = map.find(v)->first;
//        Key vMap = tree.find(v)->first;
//        assert(vMap == vRef);
//    }
    assert(tree.size() == map.size());

    return 0;
}

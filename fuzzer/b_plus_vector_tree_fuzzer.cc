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

#include <vector>

#include "include/phtree/common/bpt_vector_tree.h"

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
        std::cout << "    vector_tree<Id, 8> tree{};" << std::endl;
        std::cout << "    std::vector<Id> ref{};" << std::endl;
    }

    auto scopeguard = []() { std::cout << "};" << std::endl; };

    phtree::bptree::detail::vector_tree<Id, 8> tree;
    std::vector<Id> ref;

    size_t pos = 0;

    while (pos + 4 < Size) {
        Instruction inst = Data[pos++] % 3;
        Key key = Data[pos++];
        Value value = Data[pos++];
        switch (inst) {
        case 0: {
            if (PRINT) {
                std::cout << "    tree.emplace_back(" << (int)key << ", " << (int)value << ");"
                          << std::endl;
                std::cout << "    ref.emplace_back(" << (int)key << ", " << (int)value << ");"
                          << std::endl;
            }
            tree.emplace_back(key, value);
            ref.emplace_back(key, value);
            break;
        }
        case 1: {
            if (!tree.empty()) {
                if (PRINT) {
                    std::cout << "    tree.erase_back();" << std::endl;
                    std::cout << "    ref.erase(ref.end() - 1);" << std::endl;
                }
                tree.erase_back();
                ref.erase(ref.end() - 1);
            }
            break;
        }
        case 2: {
            if (!tree.empty()) {
                size_t index = key % tree.size();
                if (PRINT) {
                    std::cout << "    tree[" << (int)key << " % tree.size()] =  std::make_pair("
                              << (int)key << ", " << (int)value << ");" << std::endl;
                    std::cout << "    ref[" << (int)key << " % ref.size()] =  std::make_pair("
                              << (int)key << ", " << (int)value << ");" << std::endl;
                }
                tree[index] = std::make_pair(key, value);
                ref[index] = std::make_pair(key, value);
            }
            break;
        }
        default:
            std::cout << "Unexpected instruction: " << inst << std::endl;
            assert(false);
        }
    }

    //    tree._check();
    //std::cout << " sizes:   " << tree.size() << " == " << ref.size() << std::endl;
    assert(tree.size() == ref.size());

    for (size_t i = 0; i < tree.size(); ++i) {
        if (!true) {
            std::cout << "    " << (int)tree[i].first << " == " << (int)ref[i].first << std::endl;
        }
        assert(tree[i].first == ref[i].first);
        assert(tree[i].second == ref[i].second);
    }

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

    return 0;
}

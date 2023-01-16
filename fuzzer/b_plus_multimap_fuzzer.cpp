
#include <assert.h>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <ostream>

#include <map>

#include "include/phtree/common/b_plus_tree_multimap.h"

/*
 *   clang++ -g -std=c++17 -fsanitize=fuzzer include/phtree/common/BPT_MM_Fuzzer.cpp
 *   ./a.out
 *   ./a.out -artifact_prefix=/home/franky/tmp/fuzz/artifacts/
 *
 *   ./a.out -artifact_prefix=/home/franky/tmp/fuzz/artifacts/  -minimize_crash=1 -runs=10000
 *   /home/franky/tmp/fuzz/artifacts/crash-75171a275e017e3b18b3a9094203b521097a4f49
 *
 *   ./a.out -artifact_prefix=/home/franky/tmp/fuzz/artifacts/
 *   /home/franky/tmp/fuzz/artifacts/crash-75171a275e017e3b18b3a9094203b521097a4f49
 *
 *   ./a.out -artifact_prefix=/home/franky/tmp/fuzz/artifacts/
 *   /home/franky/tmp/fuzz/artifacts/minimized-from-185ecf42f208c2a7736a98ba0403f31868bcb681
 *
 *   rm -rf /home/franky/tmp/fuzz/artifacts/*
 */

static volatile int Sink;

using Instruction = std::uint8_t;
using Key = std::uint8_t;
using Value = std::uint8_t;

constexpr bool PRINT = !true;

void print() {}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* Data, size_t Size) {
    assert(Data);

    if (PRINT) {
        std::cout << "TEST(PhTreeBptMulitmapTest, FuzzTest1) {" << std::endl;
        std::cout << "    using Key = std::uint8_t;" << std::endl;
        std::cout << "    using Value = std::uint8_t;" << std::endl;
        std::cout << "    b_plus_tree_multimap<Key, Value> tree{};" << std::endl;
    }

    auto scopeguard = []() { std::cout << "};" << std::endl; };

    improbable::phtree::b_plus_tree_multimap<Key, Value> tree;
    std::multimap<Key, Value> map;

    size_t pos = 0;

    while (pos + 4 < Size) {
        Instruction inst = Data[pos++] % 4;
        Key key = Data[pos++];
        Value value = Data[pos++];
        switch (inst) {
        case 0: {
            if (PRINT)
                std::cout << "    tree.emplace(" << (int)key << ", " << (int)value << ");"
                          << std::endl;
            tree.emplace(key, value);
            map.emplace(key, value);
            break;
        }
        case 1: {
            if (PRINT)
                std::cout << "    tree.erase(" << (int)key << ");" << std::endl;
            tree.erase(key);
            map.erase(key);
            break;
        }
        case 2: {
            if (PRINT)
                std::cout << "    auto it = find.find(" << (int)key << ");" << std::endl;
            auto it = tree.find(key);
            if (it != tree.end()) {
                if (PRINT)
                    std::cout << "    tree.erase(it);" << std::endl;
                tree.erase(it);
            }
            auto it2 = map.find(key);
            if (it2 != map.end()) {
                map.erase(it2);
            }
            break;
        }
        case 3: {
            if (PRINT)
                std::cout << "    auto it = find.lower_bound(" << (int)key << ");" << std::endl;
            auto it = tree.lower_bound(key);
            if (PRINT)
                std::cout << "    tree.emplace_hint(it" << (int)key << ", " << (int)value << ");" << std::endl;
            tree.emplace_hint(it, key, value);
            auto it2 = map.lower_bound(key);
            map.emplace_hint(it2, key, value);
            break;
        }
        default:
            std::cout << "Unexpected instruction: " << inst << std::endl;
        }
    }

    tree._check();

    for (auto& entry : map) {
        const Key& vRef = entry.first;
        Key vMap = tree.find(vRef)->first;
        assert(vMap == vRef);
    }
    for (auto& entry : tree) {
        Key v = entry.first;
        const Key& vRef = map.find(v)->first;
        Key vMap = tree.find(v)->first;
        assert(vMap == vRef);
    }
    assert(tree.size() == map.size());

    return 0;
}

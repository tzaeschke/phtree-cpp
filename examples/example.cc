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

// phtree-cmake-1.cpp : Defines the entry point for the application.
//

//#include "phtree-cmake-1.h"

#include "../phtree/phtree.h"
#include "../phtree/phtree_multimap.h"
#include <chrono>
#include <iostream>
#include <set>

using namespace improbable::phtree;

template <dimension_t DIM = 2, size_t AREA_LEN = 1000, size_t LEVELS = 21>
struct ConverterWithLevels : public ConverterPointBase<DIM, double, scalar_64_t> {
    static_assert(LEVELS >= 1 && "There must be at least one level");
    static constexpr double divider_ = 1 << (LEVELS - 1);  // = 2 ^ (LEVELS - 1);
    static constexpr double multiplier_ = 1. / divider_;

    explicit ConverterWithLevels() {}

    [[nodiscard]] PhPoint<DIM> pre(const PhPointD<DIM>& point) const {
        PhPoint<DIM> out;
        for (dimension_t i = 0; i < DIM; ++i) {
            out[i] = point[i] * multiplier_;
        }
        return out;
    }

    [[nodiscard]] PhPointD<DIM> post(const PhPoint<DIM>& in) const {
        PhPointD<DIM> out;
        for (dimension_t i = 0; i < DIM; ++i) {
            out[i] = ((double)in[i]) * divider_;
        }
        return out;
    }

    [[nodiscard]] auto pre_query(const PhBoxD<DIM>& query_box) const {
        return PhBox{pre(query_box.min()), pre(query_box.max())};
    }
};

template <dimension_t DIM>
struct MyConverterMultiply : public ConverterPointBase<DIM, double, scalar_64_t> {
    explicit MyConverterMultiply(double multiplier)
    : multiplier{multiplier}, divider_{1. / multiplier} {}

    [[nodiscard]] PhPoint<DIM> pre(const PhPointD<DIM>& point) const {
        PhPoint<DIM> out;
        for (dimension_t i = 0; i < DIM; ++i) {
            out[i] = point[i] * multiplier;
        }
        return out;
    }

    [[nodiscard]] PhPointD<DIM> post(const PhPoint<DIM>& in) const {
        PhPointD<DIM> out;
        for (dimension_t i = 0; i < DIM; ++i) {
            out[i] = ((double)in[i]) * divider;
        }
        return out;
    }

    [[nodiscard]] auto pre_query(const PhBoxD<DIM>& query_box) const {
        return PhBox{pre(query_box.min()), pre(querybox.max())};
    }

    const double multiplier;
    const double divider_;
};

int main_issue_60_2() {
    //   MyConverterMultiply<2> converter{ 1000000 };
    //   auto tree = PhTreeMultiMapD<2, int, MyConverterMultiply<2>,
    //   std::unordered_set<int>>(converter);
    //auto tree = PhTreeMultiMapD<2, int, ConverterIEEE<2>, std::unordered_set<int>>();
    auto tree = PhTreeMultiMapD<2, int, ConverterWithLevels<2>, std::unordered_set<int>>();
    std::vector<PhPointD<2>> vecPos;
    int dim = 1000;

    int num = 30000;
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = {(double)(rand() % dim), (double)(rand() % dim)};
        vecPos.push_back(p);
        tree.emplace(p, i);
    }

    long T = 0;
    int nT = 0;
    while (true) {
        auto t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num; ++i) {
            PhPointD<2>& p = vecPos[i];
            PhPointD<2> newp = {(double)(rand() % dim), (double)(rand() % dim)};
            tree.relocate(p, newp, i, true);
            p = newp;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto s = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        ++nT;
        T += s.count() / 1000;
        std::cout << s.count() << "    " << (T / nT)
                  << "     msec/num= " << (s.count() / (double)num) << std::endl;
    }

    return 0;
}

int main_issue_60() {
    // auto tree = PhTreeMultiMapD<2, int>();
    PhTreeMultiMapD<2, int, ConverterIEEE<2>, std::set<int>> tree;
    std::vector<PhPointD<2>> vecPos;
    int dim = 10000;

    int num = 10000000;
    auto start1 = std::chrono::steady_clock::now();
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = {(double)(rand() % dim), (double)(rand() % dim)};
        vecPos.push_back(p);
        tree.emplace(p, i);
    }
    auto end1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
    std::cout << "elapsed time 1: " << elapsed_seconds1.count() << "s\n";

    auto start2 = std::chrono::steady_clock::now();
    for (int i = 0; i < num; ++i) {
        PhPointD<2>& p = vecPos[i];
        PhPointD<2> newp = {(double)(rand() % dim), (double)(rand() % dim)};
        tree.relocate(p, newp, i);
        p = newp;
    }
    auto end2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
    std::cout << "elapsed time 2: " << elapsed_seconds2.count() << "s\n";

    // for (int i = 0; i < num; ++i) {
    //    PhPointD<2> p = vecPos[i];
    //    for (auto it = tree.begin_knn_query(1, p, DistanceEuclidean<2>()); it != tree.end(); ++it)
    //    {
    //        std::cout << "nn: " << i << " -> " << *it << std::endl;
    //    }
    //}

    return 0;
}

int main() {
    std::cout << "PH-Tree example with 3D `double` coordinates." << std::endl;

    MyConverterMultiply<2> converter{1000000};
    auto tree = PhTreeMultiMapD<2, int, MyConverterMultiply<2>, std::unordered_set<int>>(converter);
    // PhPointD<3> p1({ 1, 1, 1 });
    // PhPointD<3> p2({ 2, 2, 2 });
    // PhPointD<3> p3({ 3, 3, 3 });
    // PhPointD<3> p4({ 4, 4, 4 });

    //// PhTreeD<3, int> tree;
    // PhTreeD<3, int> tree;
    // tree.emplace(p1, 1);
    // tree.emplace(p2, 2);
    // tree.emplace(p3, 3);
    // tree.emplace(p4, 4);

    // std::cout << "All values:" << std::endl;
    // for (auto it : tree) {
    //    std::cout << "    id=" << it << std::endl;
    //}
    // std::cout << std::endl;

    // std::cout << "All points in range:" << p2 << "/" << p4 << std::endl;
    // for (auto it = tree.begin_query({ p2, p4 }); it != tree.end(); ++it) {
    //    std::cout << "    " << it.second() << " -> " << it.first() << std::endl;
    //}
    // std::cout << std::endl;

    // std::cout << "PH-Tree is a MAP which means that, like std::map, every position " <<
    // std::endl; std::cout << " (=key) can have only ONE value." << std::endl; std::cout << "Storing
    // multiple values for a single coordinate requires storing " << std::endl; std::cout << "lists
    // or sets, for example using PhTree<3, std::vector<int>>." << std::endl;

    // PhPointD<3> p4b({ 4, 4, 4 });
    // tree.emplace(p4b, 5);
    //// Still showing '4' after emplace()
    // std::cout << "ID at " << p4b << ": " << tree.find(p4b).second() << std::endl;

    std::cout << "Done." << std::endl;

    // main_issue_60();
    main_issue_60_2();

    /*   _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
       _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
       _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
       _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
       _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
       _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
       std::cout << "leak: " << _CrtDumpMemoryLeaks() << std::endl;*/

    std::cout << "Done." << std::endl;

    return 0;
}

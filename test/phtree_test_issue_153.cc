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

#include "phtree/phtree.h"
#include "phtree/phtree_multimap.h"
#include <include/gtest/gtest.h>
#include <random>
#include <vector>
#include<iostream>
using namespace improbable::phtree;

namespace phtree_test_issue_153 {

// Number of entries that have the same coordinate
static const size_t NUM_DUPL = 1;
static const double WORLD_MIN = -100;
static const double WORLD_MAX = 100;
static const double BOX_LEN = 1;

template <dimension_t DIM>
using TestKey = PhBoxD<DIM>;

class DoubleRng {
  public:
    DoubleRng(double minIncl, double maxExcl) : eng(), rnd{minIncl, maxExcl} {}

    double next() {
        return rnd(eng);
    }

  private:
    std::default_random_engine eng;
    std::uniform_real_distribution<double> rnd;
};

template <dimension_t DIM>
void generateCube(std::vector<TestKey<DIM>>& points, size_t N) {
    assert(N % NUM_DUPL == 0);
    DoubleRng rng(WORLD_MIN, WORLD_MAX);
    auto reference_set = std::unordered_map<TestKey<DIM>, size_t>();

    points.reserve(N);
    for (size_t i = 0; i < N / NUM_DUPL; i++) {
        // create duplicates, i.e. entries with the same coordinates. However, avoid unintentional
        // duplicates.
        TestKey<DIM> key{};
        for (dimension_t d = 0; d < DIM; ++d) {
            key.min()[d] = rng.next();
            key.max()[d] = key.min()[d] + BOX_LEN;
        }
        if (reference_set.count(key) != 0) {
            i--;
            continue;
        }
        reference_set.emplace(key, i);
        for (size_t dupl = 0; dupl < NUM_DUPL; dupl++) {
            auto point = TestKey<DIM>(key);
            points.push_back(point);
        }
    }
    ASSERT_EQ(reference_set.size(), N / NUM_DUPL);
    ASSERT_EQ(points.size(), N);
}

TEST(PhTreeTestIssue153, TestIssue153) {
    /*
     * This used to cause a heap-use-after-free problem in sparse_map
     */
    auto points = std::vector<PhBoxD<2>>();
    points.emplace_back(PhBoxD<2>{{1, 1}, {1, 1}});
    points.emplace_back(PhBoxD<2>{{3, 3}, {3, 3}});
    points.emplace_back(PhBoxD<2>{{1, 3}, {1, 3}});
    points.emplace_back(PhBoxD<2>{{3, 1}, {3, 1}});
    points.emplace_back(PhBoxD<2>{{0, 0}, {0, 0}});

    // The real test is ABOVE, see "#define CALLBACK"
    auto tree = PhTreeMultiMapBoxD<2, size_t>();
    for (size_t i = 0; i < points.size(); ++i) {
        tree.emplace(points[i], i);
    }

    PhBoxD<2> b{{2, 2}, {2, 2}};
    for (size_t r = 0; r < 100; ++r) {
        for (size_t i = 0; i < points.size(); ++i) {
            PhBoxD<2>& bOld = points[i];
            double m0 = (bOld.min()[0] + b.min()[0]) / 2;
            double m1 = (bOld.min()[1] + b.min()[1]) / 2;
            double m2 = (bOld.max()[0] + b.max()[0]) / 2;
            double m3 = (bOld.max()[1] + b.max()[1]) / 2;
            PhBoxD<2> bNew{{m0, m1}, {m2, m3}};
            ASSERT_EQ(1, tree.relocate(bOld, bNew, i));
            points[i] = bNew;
        }
    }
}

/*
 * Try move in ever smaller steps.
 */
TEST(PhTreeTestIssue153, TestIssue153_2) {
    int N = 10;
    const dimension_t DIM = 2;
    std::vector<TestKey<DIM>> points;
    generateCube(points, N);

    auto tree = PhTreeMultiMapBoxD<DIM, size_t>();
    for (size_t i = 0; i < points.size(); ++i) {
        tree.emplace(points[i], i);
    }

    TestKey<DIM> x{};
    for (dimension_t d = 0; d < DIM; ++d) {
        x.min()[d] = 2.;
        x.max()[d] = 2.;
    }

    for (size_t r = 0; r < 100; ++r) {
        for (size_t i = 0; i < points.size(); ++i) {
            TestKey<DIM>& bOld = points[i];
            TestKey<DIM> bNew{};
            for (dimension_t d = 0; d < DIM; ++d) {
                double m0 = (bOld.min()[d] + x.min()[d]) / 2;
                double m1 = (bOld.max()[d] + x.max()[d]) / 2;
                bNew.min()[d] = m0;
                bNew.max()[d] = m1;
            }
            ASSERT_EQ(1, tree.relocate(bOld, bNew, i));
            points[i] = bNew;
        }
    }
}

/*
 * Try moving in very small (but const-sized) steps.
 */
TEST(PhTreeTestIssue153, TestIssue153_Linear) {
    int N = 10;
    const int R_MAX = 100000;
    const dimension_t DIM = 2;
    std::vector<TestKey<DIM>> points;
    generateCube(points, N);

    TestKey<DIM> x{};
    for (dimension_t d = 0; d < DIM; ++d) {
        x.min()[d] = 2.;
        x.max()[d] = 2. + BOX_LEN;
    }

    std::vector<PhPointD<DIM>> moves;

    auto tree = PhTreeMultiMapBoxD<DIM, size_t>();
    for (size_t i = 0; i < points.size(); ++i) {
        tree.emplace(points[i], i);
        // distances.emplace_back(distance(points[i], x));
        PhPointD<DIM>& m = moves.emplace_back();
        for (dimension_t d = 0; d < DIM; ++d) {
            m[d] = (x.min()[d] - points[i].min()[d]) / R_MAX;
        }
    }

    for (size_t r = 0; r < R_MAX; ++r) {
        for (size_t i = 0; i < points.size(); ++i) {
            TestKey<DIM>& bOld = points[i];
            if (i == 0)
                std::cout << "r=" << r << " i=" << i << " " << bOld << std::endl;
            TestKey<DIM> bNew{};
            PhPointD<DIM> mov = moves[i];
            for (dimension_t d = 0; d < DIM; ++d) {
                bNew.min()[d] = bOld.min()[d] + mov[d];
                bNew.max()[d] = bOld.max()[d] + mov[d];
            }
            ASSERT_EQ(1, tree.relocate(bOld, bNew, i));
            points[i] = bNew;
        }
    }
}

/*
 * Try moving in very small (but const-sized) steps.
 */
TEST(PhTreeTestIssue153, TestIssue153_New2025) {
 class IdentityScalarConverter {
      public:
        static std::int64_t pre(std::int64_t value) {
            return value;
        }

        static std::int64_t post(std::int64_t value) {
            return value;
        }
    };
    using multi_map_box = PhTreeMultiMapBox<
        2,
        int,
        SimpleBoxConverter<2, std::int64_t, std::int64_t, IdentityScalarConverter>>;

    multi_map_box tree;

    tree.emplace({{-4294'967296, -4294'967296}, {4294'967296, 4294'967296}}, 0);
    tree.emplace({{-4294'967296, -4294'967296}, {4294'967296, 4294'967296}}, 1);
    tree.relocate(
        {{-4294'967296, -4294'967296}, {4294'967296, 4294'967296}},
        {{-6442'450944, -4294'967296}, {2147'483648, 4294'967296}},
        1);
    tree.relocate(
        {{-6442'450944, -4294'967296}, {2147'483648, 4294'967296}},
        {{-8589'934592, -4294'967296}, {0, 4294'967296}},
        1);
    tree.relocate(
        {{-8589'934592, -4294'967296}, {0, 4294'967296}},
        {{-10737'418240, -4294'967296}, {-2147'483648, 4294'967296}},
        1);
    tree.emplace({{53687'091200, 53687'091200}, {66571'993088, 66571'993088}}, 2);
    tree.emplace({{-118111'600640, 2147'483648}, {-105226'698752, 15032'385536}}, 3);
    tree.relocate(
        {{97444'549454, 97444'549454}, {106034'484046, 106034'484046}},
        {{98963'049704, 98963'049704}, {107552'984296, 107552'984296}},
        0);
    tree.emplace({{-4294'967296, -4294'967296}, {4294'967296, 4294'967296}}, 4);
    tree.relocate(
        {{-4294'967296, -4294'967296}, {4294'967296, 4294'967296}},
        {{-4294'967296, -6442'450944}, {4294'967296, 2147'483648}},
        4);
    tree.relocate(
        {{-4294'967296, -6442'450944}, {4294'967296, 2147'483648}},
        {{-4294'967296, -8589'934592}, {4294'967296, 0}},
        4);
    PhTreeDebugHelper::CheckConsistency(tree);  // This line passes fine.
    std::cout<<"last op" << std::endl;
    tree.relocate(
        {{-4294'967296, -8589'934592}, {4294'967296, 0}},
        {{-4294'967296, -10737'418240}, {4294'967296, -2147'483648}},
        4);  // This causes a heap-use-after-free.
}

}  // namespace phtree_test_issue_153

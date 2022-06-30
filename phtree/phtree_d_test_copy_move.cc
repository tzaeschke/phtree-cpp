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

#include "phtree/phtree.h"
#include <gtest/gtest.h>
#include <random>

using namespace improbable::phtree;

template <dimension_t DIM>
using TestPoint = PhPointD<DIM>;

template <dimension_t DIM, typename T>
using TestTree = PhTreeD<DIM, T>;

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

struct IdCopyOnly {
    explicit IdCopyOnly(const size_t i) : _i{static_cast<int>(i)} {}

    IdCopyOnly() = default;
    IdCopyOnly(const IdCopyOnly& other) = default;
    IdCopyOnly(IdCopyOnly&& other) = delete;
    IdCopyOnly& operator=(const IdCopyOnly& other) = default;
    IdCopyOnly& operator=(IdCopyOnly&& other) = delete;
    ~IdCopyOnly() = default;

    bool operator==(const IdCopyOnly& rhs) const {
        return _i == rhs._i;
    }

    int _i;
};

struct IdMoveOnly {
    explicit IdMoveOnly(const size_t i) : _i{static_cast<int>(i)} {}

    IdMoveOnly() = default;
    IdMoveOnly(const IdMoveOnly& other) = delete;
    IdMoveOnly(IdMoveOnly&& other) = default;
    IdMoveOnly& operator=(const IdMoveOnly& other) = delete;
    IdMoveOnly& operator=(IdMoveOnly&& other) = default;
    ~IdMoveOnly() = default;

    bool operator==(const IdMoveOnly& rhs) const {
        return _i == rhs._i;
    }

    int _i;
};

// Assert that copy-ctr is not called even when available
struct IdCopyOrMove {
    explicit IdCopyOrMove(const size_t i) : _i{static_cast<int>(i)} {}

    IdCopyOrMove() = default;
    IdCopyOrMove(const IdCopyOrMove&) {
        assert(false);
    }
    IdCopyOrMove(IdCopyOrMove&& other) = default;
    IdCopyOrMove& operator=(const IdCopyOrMove&) {
        assert(false);
    }
    IdCopyOrMove& operator=(IdCopyOrMove&& other) = default;
    ~IdCopyOrMove() = default;

    bool operator==(const IdMoveOnly& rhs) const {
        return _i == rhs._i;
    }

    int _i;
};

struct PointDistance {
    PointDistance(double distance, size_t id) : _distance(distance), _id(id) {}

    double _distance;
    size_t _id;
};

bool comparePointDistance(PointDistance& i1, PointDistance& i2) {
    return (i1._distance < i2._distance);
}

template <dimension_t DIM>
double distance(const TestPoint<DIM>& p1, const TestPoint<DIM>& p2) {
    double sum2 = 0;
    for (dimension_t i = 0; i < DIM; i++) {
        double d = p1[i] - p2[i];
        sum2 += d * d;
    }
    return sqrt(sum2);
}

template <dimension_t DIM>
double distanceL1(const TestPoint<DIM>& p1, const TestPoint<DIM>& p2) {
    double sum = 0;
    for (dimension_t i = 0; i < DIM; i++) {
        sum += std::abs(p1[i] - p2[i]);
    }
    return sum;
}

template <dimension_t DIM>
void generateCube(std::vector<TestPoint<DIM>>& points, size_t N) {
    DoubleRng rng(-1000, 1000);
    auto refTree = std::map<TestPoint<DIM>, size_t>();

    points.reserve(N);
    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM> point{};
        for (dimension_t d = 0; d < DIM; ++d) {
            point[d] = rng.next();
        }
        if (refTree.count(point) != 0) {
            i--;
            continue;
        }

        refTree.emplace(point, i);
        points.push_back(point);
    }
    ASSERT_EQ(refTree.size(), N);
    ASSERT_EQ(points.size(), N);
}

template <dimension_t DIM, typename Id>
void SmokeTestBasicOps(size_t N) {
    TestTree<DIM, Id> tree;
    std::vector<TestPoint<DIM>> points;
    generateCube(points, N);

    ASSERT_EQ(0, tree.size());
    ASSERT_TRUE(tree.empty());
    PhTreeDebugHelper::CheckConsistency(tree);

    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM>& p = points.at(i);
        ASSERT_EQ(tree.count(p), 0);
        ASSERT_EQ(tree.end(), tree.find(p));

        Id id(i);
        //        if (i % 2 == 0) {
        // ASSERT_TRUE(tree.try_emplace(p, id).second);
        //        } else {
        ASSERT_TRUE(tree.insert(p, id).second);
        //        }
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_NE(tree.end(), tree.find(p));
        ASSERT_EQ(id._i, tree.find(p)->_i);
        ASSERT_EQ(i + 1, tree.size());

        // try adding it again
        ASSERT_FALSE(tree.insert(p, id).second);
        ASSERT_FALSE(tree.emplace(p, id).second);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_NE(tree.end(), tree.find(p));
        ASSERT_EQ(id._i, tree.find(p)->_i);
        ASSERT_EQ(i + 1, tree.size());
        ASSERT_FALSE(tree.empty());
    }

    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM>& p = points.at(i);
        auto q = tree.begin_query({p, p});
        ASSERT_NE(q, tree.end());
        ASSERT_EQ(i, (*q)._i);
        q++;
        ASSERT_EQ(q, tree.end());
    }

    PhTreeDebugHelper::CheckConsistency(tree);

    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM>& p = points.at(i);
        ASSERT_NE(tree.find(p), tree.end());
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(i, tree.find(p)->_i);
        if (i % 2 == 0) {
            ASSERT_EQ(1, tree.erase(p));
        } else {
            auto iter = tree.find(p);
            ASSERT_EQ(1, tree.erase(iter));
        }

        ASSERT_EQ(tree.count(p), 0);
        ASSERT_EQ(tree.end(), tree.find(p));
        ASSERT_EQ(N - i - 1, tree.size());

        // try remove again
        ASSERT_EQ(0, tree.erase(p));
        ASSERT_EQ(tree.count(p), 0);
        ASSERT_EQ(tree.end(), tree.find(p));
        ASSERT_EQ(N - i - 1, tree.size());
        if (i < N - 1) {
            ASSERT_FALSE(tree.empty());
        }
    }
    ASSERT_EQ(0, tree.size());
    ASSERT_TRUE(tree.empty());
    PhTreeDebugHelper::CheckConsistency(tree);
}

TEST(PhTreeDTestCopyMove, SmokeTestBasicOpsCopyOnly) {
    SmokeTestBasicOps<1, IdCopyOnly>(100);
    SmokeTestBasicOps<3, IdCopyOnly>(100);
    SmokeTestBasicOps<6, IdCopyOnly>(100);
    SmokeTestBasicOps<10, IdCopyOnly>(100);
    SmokeTestBasicOps<20, IdCopyOnly>(100);
    SmokeTestBasicOps<63, IdCopyOnly>(100);
}

template <dimension_t DIM, typename Id>
void SmokeTestBasicOpsMoveOnly(size_t N) {
    TestTree<DIM, Id> tree;
    std::vector<TestPoint<DIM>> points;
    generateCube(points, N);

    ASSERT_EQ(0, tree.size());
    ASSERT_TRUE(tree.empty());
    PhTreeDebugHelper::CheckConsistency(tree);

    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM>& p = points.at(i);
        ASSERT_EQ(tree.count(p), 0);
        ASSERT_EQ(tree.end(), tree.find(p));

        if (i % 2 == 0) {
            ASSERT_TRUE(tree.try_emplace(p, Id(i)).second);
        } else {
            ASSERT_TRUE(tree.emplace(p, Id(i)).second);
        }
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_NE(tree.end(), tree.find(p));
        ASSERT_EQ(i, tree.find(p)->_i);
        ASSERT_EQ(i + 1, tree.size());

        // try adding it again
        ASSERT_FALSE(tree.emplace(p, Id(i)).second);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_NE(tree.end(), tree.find(p));
        ASSERT_EQ(i, tree.find(p)->_i);
        ASSERT_EQ(i + 1, tree.size());
        ASSERT_FALSE(tree.empty());
    }

    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM>& p = points.at(i);
        auto q = tree.begin_query({p, p});
        ASSERT_NE(q, tree.end());
        ASSERT_EQ(i, (*q)._i);
        q++;
        ASSERT_EQ(q, tree.end());
    }

    PhTreeDebugHelper::CheckConsistency(tree);

    for (size_t i = 0; i < N; i++) {
        TestPoint<DIM>& p = points.at(i);
        ASSERT_NE(tree.find(p), tree.end());
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(i, tree.find(p)->_i);
        if (i % 2 == 0) {
            ASSERT_EQ(1, tree.erase(p));
        } else {
            auto iter = tree.find(p);
            ASSERT_EQ(1, tree.erase(iter));
        }

        ASSERT_EQ(tree.count(p), 0);
        ASSERT_EQ(tree.end(), tree.find(p));
        ASSERT_EQ(N - i - 1, tree.size());

        // try remove again
        ASSERT_EQ(0, tree.erase(p));
        ASSERT_EQ(tree.count(p), 0);
        ASSERT_EQ(tree.end(), tree.find(p));
        ASSERT_EQ(N - i - 1, tree.size());
        if (i < N - 1) {
            ASSERT_FALSE(tree.empty());
        }
    }
    ASSERT_EQ(0, tree.size());
    ASSERT_TRUE(tree.empty());
    PhTreeDebugHelper::CheckConsistency(tree);
}

TEST(PhTreeDTestCopyMove, SmokeTestBasicOpsMoveOnly) {
    SmokeTestBasicOpsMoveOnly<1, IdMoveOnly>(100);
    SmokeTestBasicOpsMoveOnly<3, IdMoveOnly>(100);
    SmokeTestBasicOpsMoveOnly<6, IdMoveOnly>(100);
    SmokeTestBasicOpsMoveOnly<10, IdMoveOnly>(100);
    SmokeTestBasicOpsMoveOnly<20, IdMoveOnly>(100);
    SmokeTestBasicOpsMoveOnly<63, IdMoveOnly>(100);
}

TEST(PhTreeDTestCopyMove, SmokeTestBasicOpsCopyFails) {
    SmokeTestBasicOpsMoveOnly<1, IdCopyOrMove>(100);
    SmokeTestBasicOpsMoveOnly<3, IdCopyOrMove>(100);
    SmokeTestBasicOpsMoveOnly<6, IdCopyOrMove>(100);
    SmokeTestBasicOpsMoveOnly<10, IdCopyOrMove>(100);
    SmokeTestBasicOpsMoveOnly<20, IdCopyOrMove>(100);
    SmokeTestBasicOpsMoveOnly<63, IdCopyOrMove>(100);
}

TEST(PhTreeDTest, TestInsert) {
    using Id = IdCopyOnly;
    const dimension_t dim = 3;
    TestTree<dim, Id> tree;
    size_t N = 1000;

    std::vector<TestPoint<dim>> points;
    generateCube(points, N);

    for (size_t i = 0; i < N; i++) {
        TestPoint<dim>& p = points.at(i);
        Id id(i);
        ASSERT_EQ(true, tree.insert(p, id).second);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(id._i, tree.find(p)->_i);

        // try add again
        ASSERT_EQ(false, tree.insert(p, id).second);
        ASSERT_EQ(i, tree.insert(p, id).first._i);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(id._i, tree.find(p)->_i);
    }
    ASSERT_EQ(N, tree.size());

    for (size_t i = 0; i < N; i++) {
        TestPoint<dim>& p = points.at(i);
        auto q = tree.begin_query({p, p});
        ASSERT_NE(q, tree.end());
        ASSERT_EQ(i, (*q)._i);
        q++;
        ASSERT_EQ(q, tree.end());
    }

    for (size_t i = 0; i < N; i++) {
        TestPoint<dim>& p = points.at(i);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(i, tree.find(p)->_i);
    }
}

TEST(PhTreeDTest, TestEmplace) {
    using Id = IdMoveOnly;
    const dimension_t dim = 3;
    TestTree<dim, Id> tree;
    size_t N = 1000;

    std::vector<TestPoint<dim>> points;
    generateCube(points, N);

    for (size_t i = 0; i < N; i++) {
        TestPoint<dim>& p = points.at(i);
        ASSERT_EQ(true, tree.emplace(p, Id(i)).second);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(i, tree.find(p)->_i);
        ASSERT_EQ(i + 1, tree.size());

        // try adding again, this should _not_ replace the existing value
        ASSERT_EQ(false, tree.emplace(p, Id(-i)).second);
        ASSERT_EQ(i, tree.emplace(p, Id(i)).first._i);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(i, tree.find(p)->_i);

        // Check that the returned value is a reference
        tree.emplace(p, Id(-i)).first._i++;
        ASSERT_EQ(i + 1, tree.emplace(p, Id(i)).first._i);
        tree.emplace(p, Id(-i)).first = Id(i);
        ASSERT_EQ(i, tree.emplace(p, Id(i)).first._i);
    }
    ASSERT_EQ(N, tree.size());

    for (size_t i = 0; i < N; i++) {
        TestPoint<dim>& p = points.at(i);
        auto q = tree.begin_query({p, p});
        ASSERT_NE(q, tree.end());
        ASSERT_EQ(i, (*q)._i);
        q++;
        ASSERT_EQ(q, tree.end());
    }

    for (size_t i = 0; i < N; i++) {
        TestPoint<dim>& p = points.at(i);
        ASSERT_EQ(tree.count(p), 1);
        ASSERT_EQ(i, tree.find(p)->_i);
    }
}

// TEST(PhTreeDTest, TestSquareBrackets) {
//     const dimension_t dim = 3;
//     TestTree<dim, Id> tree;
//     size_t N = 1000;
//
//     std::vector<TestPoint<dim>> points;
//     generateCube(points, N);
//
//     for (size_t i = 0; i < N; i++) {
//         TestPoint<dim>& p = points.at(i);
//         Id id(i);
//         ASSERT_EQ(0, tree[p]._i);
//         ASSERT_EQ(tree.count(p), 1);
//         if (i % 2 == 0) {
//             tree[p]._i = i;
//         } else {
//             tree[p] = id;
//         }
//         ASSERT_EQ(id._i, tree.find(p)->_i);
//         ASSERT_EQ(i + 1, tree.size());
//
//         // try `add` again
//         ASSERT_EQ(i, tree[p]._i);
//         ASSERT_EQ(tree.count(p), 1);
//         ASSERT_EQ(id._i, tree.find(p)->_i);
//     }
//     ASSERT_EQ(N, tree.size());
//
//     for (size_t i = 0; i < N; i++) {
//         TestPoint<dim>& p = points.at(i);
//         auto q = tree.begin_query({p, p});
//         ASSERT_NE(q, tree.end());
//         ASSERT_EQ(i, (*q)._i);
//         q++;
//         ASSERT_EQ(q, tree.end());
//     }
//
//     for (size_t i = 0; i < N; i++) {
//         TestPoint<dim>& p = points.at(i);
//         ASSERT_EQ(tree.count(p), 1);
//         ASSERT_EQ(i, tree.find(p)->_i);
//         ASSERT_EQ(i, tree[p]._i);
//     }
// }

// TODO relocate()
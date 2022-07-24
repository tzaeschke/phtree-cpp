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
#include "phtree_multimap.h"
#include <gtest/gtest.h>
#include <iostream>
#include <chrono>

using namespace improbable::phtree;

TEST(PhTreeTestIssues, TestIssue60) {
    auto tree = PhTreeMultiMapD<2, int>();
    std::vector<PhPointD<2>> vecPos;
    int dim = 1000;

    int num = 10;
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = { (double)(rand() % dim), (double)(rand() % dim) };
        vecPos.push_back(p);
        tree.emplace(p, i);
    }

    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = vecPos[i];
        PhPointD<2> newp = { (double)(rand() % dim), (double)(rand() % dim) };
        tree.relocate(p, newp, i);
    }
}

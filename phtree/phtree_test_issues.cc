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
#include "phtree.h"
#include "phtree_multimap.h"
#include <gtest/gtest.h>
#include <iostream>
#include <chrono>
#include <fstream>

using namespace improbable::phtree;


using namespace std;

#if defined(__clang__) || defined(__GNUC__)

void mem_usage(double &vm_usage, double &resident_set) {
    vm_usage = 0.0;
    resident_set = 0.0;
    ifstream stat_stream("/proc/self/stat", ios_base::in); //get info from proc directory
    //create some variables to get info
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;
    unsigned long vsize;
    long rss;
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
    stat_stream.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

void print_mem() {
    double vm, rss;
    mem_usage(vm, rss);
    cout << "  Virtual Memory: " << vm << " KB" << std::endl << "  Resident set size: " << rss << " KB" << endl;
}

#elif defined(_MSC_VER)
void print_mem() {
    double vm, rss;
    //mem_usage(vm, rss);
    cout << "  Virtual Memory: " << vm << " KB" << std::endl << "  Resident set size: " << rss << " KB" << endl;
}
#endif


TEST(PhTreeTestIssues, TestIssue60) {
    auto tree = PhTreeMultiMapD<2, int>();
    //auto tree = PhTreeMultiMapD<2, int, ConverterIEEE<2>, std::set<int>>();
    std::vector<PhPointD<2>> vecPos;
    int dim = 10000;
    int num = 1000000;

    auto start1 = std::chrono::steady_clock::now();
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = {(double) (rand() % dim), (double) (rand() % dim)};
        vecPos.push_back(p);
        tree.emplace(p, i);
    }

    auto end1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
    std::cout << "elapsed time 1: " << elapsed_seconds1.count() << "s\n";

    print_mem();
    auto start2 = std::chrono::steady_clock::now();
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = vecPos[i];
        PhPointD<2> newp = {(double) (rand() % dim), (double) (rand() % dim)};
        tree.relocate(p, newp, i);
    }

    auto end2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
    std::cout << "elapsed time 2: " << elapsed_seconds2.count() << "s\n";


    print_mem();
}


TEST(PhTreeTestIssues, TestIssue60_3) {
    //auto tree = PhTreeMultiMapD<3, int>();
    auto tree = PhTreeMultiMapD<3, int, ConverterIEEE<3>, std::set<int>>();
    std::vector<PhPointD<3>> vecPos;
    int dim = 10000;

    int num = 100000;
    for (int i = 0; i < num; ++i) {
        PhPointD<3> p = {(double) (rand() % dim), (double) (rand() % dim), (double) (rand() % dim)};
        vecPos.push_back(p);
        tree.emplace(p, i);
    }

    print_mem();
    for (int i = 0; i < num; ++i) {
        PhPointD<3> p = vecPos[i];
        PhPointD<3> newp = {(double) (rand() % dim), (double) (rand() % dim), (double) (rand() % dim)};
        tree.relocate(p, newp, i);
    }
    print_mem();
    FAIL();
}

TEST(PhTreeTestIssues, TestIssue6_3_MAP) {
    auto tree = PhTreeD<2, int>();
    std::vector<PhPointD<2>> vecPos;
    int dim = 10000;

    int num = 100000;
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = {(double) (rand() % dim), (double) (rand() % dim)};
        vecPos.push_back(p);
        tree.emplace(p, i);
    }

    print_mem();
    for (int i = 0; i < num; ++i) {
        PhPointD<2> p = vecPos[i];
        PhPointD<2> newp = {(double) (rand() % dim), (double) (rand() % dim)};
        tree.relocate(p, newp);
    }
    print_mem();
}



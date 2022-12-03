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
#include "benchmark_util.h"
#include "logging.h"
#include "phtree/phtree.h"
#include <benchmark/benchmark.h>
#include <random>

using namespace improbable;
using namespace improbable::phtree;
using namespace improbable::phtree::phbenchmark;

namespace {

const int GLOBAL_MAX = 10000;

enum QueryType {
    FIND,
    COUNT,
};

/*
 * Benchmark for looking up entries by their key.
 */
template <dimension_t DIM>
class IndexBenchmark {
  public:
    IndexBenchmark(
        benchmark::State& state, TestGenerator data_type, int num_entities, QueryType query_type);

    void Benchmark(benchmark::State& state);

  private:
    void SetupWorld(benchmark::State& state);
    int QueryWorldCount(benchmark::State& state);
    int QueryWorldFind(benchmark::State& state);

    const TestGenerator data_type_;
    const size_t num_entities_;
    const QueryType query_type_;

    PhTree<DIM, int> tree_;
    std::default_random_engine random_engine_;
    std::uniform_int_distribution<> cube_distribution_;
    std::vector<PhPoint<DIM>> points_;
};

template <dimension_t DIM>
IndexBenchmark<DIM>::IndexBenchmark(
    benchmark::State& state, TestGenerator data_type, int num_entities, QueryType query_type)
: data_type_{data_type}
, num_entities_(num_entities)
, query_type_(query_type)
, random_engine_{1}
, cube_distribution_{0, GLOBAL_MAX}
, points_(num_entities) {
    logging::SetupDefaultLogging();
    SetupWorld(state);
}

template <dimension_t DIM>
void IndexBenchmark<DIM>::Benchmark(benchmark::State& state) {
    int num_inner = 0;
    int num_found = 0;
    switch (query_type_) {
    case COUNT: {
        for (auto _ : state) {
            num_found += QueryWorldCount(state);
            ++num_inner;
        }
        break;
    }
    case FIND: {
        for (auto _ : state) {
            num_found += QueryWorldFind(state);
            ++num_inner;
        }
        break;
    }
    }
    // Moved outside of the loop because EXPENSIVE
    state.counters["total_result_count"] += num_found;
    state.counters["query_rate"] += num_inner;
    state.counters["result_rate"] += num_found;
    state.counters["avg_result_count"] += num_found;
}

template <dimension_t DIM>
void IndexBenchmark<DIM>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    CreatePointData<DIM>(points_, data_type_, num_entities_, 0, GLOBAL_MAX);
    for (size_t i = 0; i < num_entities_; ++i) {
        tree_.emplace(points_[i], (int)i);
    }

    state.counters["total_result_count"] = benchmark::Counter(0);
    state.counters["query_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["result_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["avg_result_count"] = benchmark::Counter(0, benchmark::Counter::kAvgIterations);

    logging::info("World setup complete.");
}

template <dimension_t DIM>
int IndexBenchmark<DIM>::QueryWorldCount(benchmark::State&) {
    static int pos = 0;
    pos = (pos + 1) % num_entities_;
    bool found = true;
    if (pos % 2 == 0) {
        assert(tree_.find(points_.at(pos)) != tree_.end());
    } else {
        int x = pos % GLOBAL_MAX;
        PhPoint<DIM> p = PhPoint<DIM>({x, x, x});
        found = tree_.find(p) != tree_.end();
    }
    return found;
}

template <dimension_t DIM>
int IndexBenchmark<DIM>::QueryWorldFind(benchmark::State&) {
    static int pos = 0;
    pos = (pos + 1) % num_entities_;
    bool found;
    if (pos % 2 == 0) {
        // This should always be a match
        found = tree_.find(points_.at(pos)) != tree_.end();
        assert(found);
    } else {
        // This should rarely be a match
        int x = pos % GLOBAL_MAX;
        PhPoint<DIM> p = PhPoint<DIM>({x, x, x});
        found = tree_.find(p) != tree_.end();
    }
    return found;
}

}  // namespace

template <typename... Arguments>
void PhTree3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_generator, num_entities, function_to_call
// PhTree 3D CUBE
BENCHMARK_CAPTURE(PhTree3D, COUNT_CU_1K, TestGenerator::CUBE, 1000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, COUNT_CU_10K, TestGenerator::CUBE, 10000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, COUNT_CU_100K, TestGenerator::CUBE, 100000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, COUNT_CU_1M, TestGenerator::CUBE, 1000000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CU_1K, TestGenerator::CUBE, 1000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CU_10K, TestGenerator::CUBE, 10000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CU_100K, TestGenerator::CUBE, 100000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CU_1M, TestGenerator::CUBE, 1000000, FIND)
    ->Unit(benchmark::kMillisecond);

// index type, scenario name, data_generator, num_entities, function_to_call
// PhTree 3D CLUSTER
BENCHMARK_CAPTURE(PhTree3D, COUNT_CL_1K, TestGenerator::CLUSTER, 1000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, COUNT_CL_10K, TestGenerator::CLUSTER, 10000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, COUNT_CL_100K, TestGenerator::CLUSTER, 100000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, COUNT_CL_1M, TestGenerator::CLUSTER, 1000000, COUNT)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CL_1K, TestGenerator::CLUSTER, 1000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CL_10K, TestGenerator::CLUSTER, 10000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CL_100K, TestGenerator::CLUSTER, 100000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, FIND_CL_1M, TestGenerator::CLUSTER, 1000000, FIND)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

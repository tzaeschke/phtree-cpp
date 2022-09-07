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

/*
 * Benchmark for window queries.
 */
template <dimension_t DIM>
class IndexBenchmark {
  public:
    IndexBenchmark(
        benchmark::State& state,
        TestGenerator data_type,
        int num_entities,
        double avg_query_result_size_);

    void Benchmark(benchmark::State& state);

  private:
    void SetupWorld(benchmark::State& state);

    void QueryWorld(benchmark::State& state, PhBox<DIM>& query);

    void CreateQuery(PhBox<DIM>& query);

    const TestGenerator data_type_;
    const size_t num_entities_;
    const double avg_query_result_size_;

    constexpr int query_endge_length() {
        return (int)(GLOBAL_MAX * pow(avg_query_result_size_ / (double)num_entities_, 1. / (double)DIM));
    };

    PhTree<DIM, int> tree_;
    std::default_random_engine random_engine_;
    std::uniform_int_distribution<> cube_distribution_;
    std::vector<PhPoint<DIM>> points_;
};

template <dimension_t DIM>
IndexBenchmark<DIM>::IndexBenchmark(
    benchmark::State& state,
    TestGenerator data_type,
    int num_entities,
    double avg_query_result_size)
: data_type_{data_type}
, num_entities_(num_entities)
, avg_query_result_size_(avg_query_result_size)
, random_engine_{1}
, cube_distribution_{0, GLOBAL_MAX}
, points_(num_entities) {
    logging::SetupDefaultLogging();
    SetupWorld(state);
}

template <dimension_t DIM>
void IndexBenchmark<DIM>::Benchmark(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        PhBox<DIM> query_box;
        CreateQuery(query_box);
        state.ResumeTiming();

        QueryWorld(state, query_box);
    }
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
void IndexBenchmark<DIM>::QueryWorld(benchmark::State& state, PhBox<DIM>& query_box) {
    int n = 0;
    for (auto q = tree_.begin_query(query_box); q != tree_.end(); ++q) {
        ++n;
    }

    state.counters["total_result_count"] += n;
    state.counters["query_rate"] += 1;
    state.counters["result_rate"] += n;
    state.counters["avg_result_count"] += n;
}

template <dimension_t DIM>
void IndexBenchmark<DIM>::CreateQuery(PhBox<DIM>& query_box) {
    int length = query_endge_length();
    // scale to ensure query lies within boundary
    double scale = (GLOBAL_MAX - (double)length) / GLOBAL_MAX;
    for (dimension_t d = 0; d < DIM; ++d) {
        scalar_64_t s = cube_distribution_(random_engine_);
        s = (scalar_64_t)(s * scale);
        query_box.min()[d] = s;
        query_box.max()[d] = s + length;
    }
}

}  // namespace

template <typename... Arguments>
void PhTree3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_type, num_entities, query_result_size
// PhTree 3D CUBE
BENCHMARK_CAPTURE(PhTree3D, WQ_CU_100_of_1K, TestGenerator::CUBE, 1000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, WQ_CU_100_of_10K, TestGenerator::CUBE, 10000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, WQ_CU_100_of_100K, TestGenerator::CUBE, 100000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, WQ_CU_100_of_1M, TestGenerator::CUBE, 1000000, 100.0)
    ->Unit(benchmark::kMillisecond);

// index type, scenario name, data_type, num_entities, query_result_size
// PhTree 3D CLUSTER
BENCHMARK_CAPTURE(PhTree3D, WQ_CL_100_of_1K, TestGenerator::CLUSTER, 1000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, WQ_CL_100_of_10K, TestGenerator::CLUSTER, 10000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, WQ_CL_100_of_100K, TestGenerator::CLUSTER, 100000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, WQ_CL_100_of_1M, TestGenerator::CLUSTER, 1000000, 100.0)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

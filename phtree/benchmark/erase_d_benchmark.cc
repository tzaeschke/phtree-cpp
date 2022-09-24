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
#include "logging.h"
#include "phtree/benchmark/benchmark_util.h"
#include "phtree/phtree.h"
#include <benchmark/benchmark.h>
#include <random>

using namespace improbable;
using namespace improbable::phtree;
using namespace improbable::phtree::phbenchmark;

namespace {

const double GLOBAL_MAX = 10000;

enum OpType { ERASE_BY_KEY, ERASE_BY_ITER };

template <dimension_t DIM>
using PointType = PhPointD<DIM>;

template <dimension_t DIM>
using TreeType = PhTreeD<DIM, size_t>;

/*
 * Benchmark for removing entries.
 */
template <dimension_t DIM, OpType OP_TYPE>
class IndexBenchmark {
  public:
    IndexBenchmark(benchmark::State& state);

    void Benchmark(benchmark::State& state);

  private:
    void SetupWorld(benchmark::State& state);
    void Insert(benchmark::State& state);
    void Remove(benchmark::State& state);

    const TestGenerator data_type_;
    const size_t num_entities_;

    TreeType<DIM> tree_;
    std::default_random_engine random_engine_;
    std::uniform_real_distribution<> cube_distribution_;
    std::vector<PhPointD<DIM>> points_;
};

template <dimension_t DIM, OpType OP_TYPE>
IndexBenchmark<DIM, OP_TYPE>::IndexBenchmark(benchmark::State& state)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, random_engine_{1}
, cube_distribution_{0, GLOBAL_MAX}
, points_(num_entities_) {
    logging::SetupDefaultLogging();
    SetupWorld(state);
}

template <dimension_t DIM, OpType OP_TYPE>
void IndexBenchmark<DIM, OP_TYPE>::Benchmark(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        auto* tree = new PhTreeD<DIM, int>();
        Insert(state);
        state.ResumeTiming();

        Remove(state);

        state.PauseTiming();
        // avoid measuring deallocation
        delete tree;
        state.ResumeTiming();
    }
}

template <dimension_t DIM, OpType OP_TYPE>
void IndexBenchmark<DIM, OP_TYPE>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    CreatePointData<DIM>(points_, data_type_, num_entities_, 0, GLOBAL_MAX);

    state.counters["total_remove_count"] = benchmark::Counter(0);
    state.counters["remove_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    logging::info("World setup complete.");
}

template <dimension_t DIM, OpType OP_TYPE>
void IndexBenchmark<DIM, OP_TYPE>::Insert(benchmark::State&) {
    for (size_t i = 0; i < num_entities_; ++i) {
        tree_.emplace(points_[i], i);
    }
}

template <dimension_t DIM>
size_t EraseByKey(TreeType<DIM>& tree, std::vector<PhPointD<DIM>>& points) {
    size_t n = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        n += tree.erase(points[i]);
    }
    return n;
}

template <dimension_t DIM>
size_t EraseByIter(TreeType<DIM>& tree, std::vector<PhPointD<DIM>>& points) {
    size_t n = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        auto iter = tree.find(points[i]);
        n += tree.erase(iter);
    }
    return n;
}

template <dimension_t DIM, OpType OP_TYPE>
void IndexBenchmark<DIM, OP_TYPE>::Remove(benchmark::State& state) {
    size_t n = 0;
    switch (OP_TYPE) {
    case OpType::ERASE_BY_KEY:
        n = EraseByKey(tree_, points_);
        break;
    case OpType::ERASE_BY_ITER:
        n = EraseByIter(tree_, points_);
        break;
    }

    state.counters["total_remove_count"] += n;
    state.counters["remove_rate"] += n;
}

}  // namespace

template <typename... Arguments>
void PhTreeEraseKey3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, OpType::ERASE_BY_KEY> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeEraseIter3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, OpType::ERASE_BY_ITER> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_type, num_entities
// PhTree with erase()
BENCHMARK(PhTreeEraseKey3D)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

// PhTree with erase(iter)
BENCHMARK(PhTreeEraseIter3D)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

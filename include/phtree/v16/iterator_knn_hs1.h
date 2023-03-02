/*
 * Copyright 2020 Improbable Worlds Limited
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

#ifndef PHTREE_V16_QUERY_KNN_HS1_H
#define PHTREE_V16_QUERY_KNN_HS1_H

#include "iterator_base.h"
#include "phtree/common/bpt_priority_queue.h"
#include "phtree/common/common.h"
#include <iostream>
#include <queue>

namespace improbable::phtree::v16 {

/*
 * kNN query implementation that uses preprocessors and distance functions.
 *
 * Implementation after Hjaltason and Samet (with some deviations: no MinDist or MaxDist used).
 * G. R. Hjaltason and H. Samet., "Distance browsing in spatial databases.", ACM TODS
 * 24(2):265--318. 1999
 */

namespace {
template <dimension_t DIM, typename T, typename SCALAR>
using EntryDist1 = std::pair<double, const Entry<DIM, T, SCALAR>*>;

template <typename ENTRY>
struct CompareEntryDistByDistance1 {
    bool operator()(const ENTRY& left, const ENTRY& right) const {
        return left.first > right.first;
    };
};
}  // namespace

template <typename T, typename CONVERT, typename DISTANCE, typename FILTER>
class IteratorKnnHS1 : public IteratorWithFilter<T, CONVERT, FILTER> {
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using KeyExternal = typename CONVERT::KeyExternal;
    using KeyInternal = typename CONVERT::KeyInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = typename IteratorWithFilter<T, CONVERT, FILTER>::EntryT;
    using EntryDistT = EntryDist1<DIM, T, SCALAR>;

  public:
    template <typename DIST, typename F>
    explicit IteratorKnnHS1(
        const EntryT& root,
        size_t min_results,
        const KeyInternal& center,
        const CONVERT* converter,
        DIST&& dist,
        F&& filter)
    : IteratorWithFilter<T, CONVERT, F>(converter, std::forward<F>(filter))
    , center_{center}
    , center_post_{converter->post(center)}
    , current_distance_{std::numeric_limits<double>::max()}
    , num_found_results_(0)
    , num_requested_results_(min_results)
    , distance_(std::forward<DIST>(dist)) {
        if (min_results <= 0 || root.GetNode().GetEntryCount() == 0) {
            this->SetFinished();
            return;
        }

        // Initialize queue, use d=0 because every imaginable point lies inside the root Node
        assert(root.IsNode());
        queue_n_.emplace(EntryDistT{0, &const_cast<EntryT&>(root)});  // TODO remove const_casts etc

        FindNextElement();
        ++N_CREATE;
        if (N_CREATE % 10000 == 0) {
            std::cout << "KNN1: " << MAX_DEPTH << " N_Q=" << N_CREATE
                      << " N_PR=" << N_PROCESSED / N_CREATE << " N_PR_R=" << N_PR_RESULT / N_CREATE
                      << " N_PR_N=" << N_PR_NODES / N_CREATE << " N_Q_R=" << N_Q_RESULT / N_CREATE
                      << " N_Q_N=" << N_Q_NODES / N_CREATE << " N_Q_N=" << N_Q_NODES_0 / N_CREATE
                      << " avg_D=" << (TOTAL_DEPTH / N_CREATE) << std::endl;
            MAX_DEPTH = 0;
            N_CREATE = 0;
            N_PROCESSED = 0;
            N_PR_RESULT = 0;
            N_PR_NODES = 0;
            N_Q_RESULT = 0;
            N_Q_NODES = 0;
            N_Q_NODES_0 = 0;
            TOTAL_DEPTH = 0;
        }
    }

    inline static size_t MAX_DEPTH{0};
    inline static long N_CREATE{0};
    inline static long N_PROCESSED{0};
    inline static long N_PR_RESULT{0};
    inline static long N_PR_NODES{0};
    inline static long N_Q_RESULT{0};
    inline static long N_Q_NODES{0};
    inline static long N_Q_NODES_0{0};
    inline static double TOTAL_DEPTH{0};

    [[nodiscard]] double distance() const {
        return current_distance_;
    }

    IteratorKnnHS1& operator++() noexcept {
        FindNextElement();
        return *this;
    }

    [[deprecated]]  // This iterator is MUCH slower!
    IteratorKnnHS1
    operator++(int) noexcept {
        IteratorKnnHS1 iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void FindNextElement() {
        while (num_found_results_ < num_requested_results_ &&
               !(queue_n_.empty() && queue_v_.empty())) {
            bool use_v = !queue_v_.empty();
            if (use_v && !queue_n_.empty()) {
                use_v = queue_n_.top() >= queue_v_.top();
            }
            ++N_PROCESSED;
            if (use_v) {
                ++N_PR_RESULT;
                // data entry
                auto& cand_v = queue_v_.top();
                ++num_found_results_;
                this->SetCurrentResult(cand_v.second);
                current_distance_ = cand_v.first;
                // We need to pop() AFTER we processed the value, otherwise the reference is
                // overwritten.
                queue_v_.pop();
                return;
            } else {
                ++N_PR_NODES;
                // inner node
                auto top = queue_n_.top();
                auto& node = top.second->GetNode();
                auto d_node = top.first;
                queue_n_.pop();

                if (queue_v_.size() >= num_requested_results_ && d_node > max_node_dist_) {
                    // ignore this node
                    continue;
                }

                for (auto& entry : node.Entries()) {
                    // auto& e2 = const_cast<EntryT&>(entry.second);
                    const auto& e2 = entry.second;
                    if (this->ApplyFilter(e2)) {
                        if (e2.IsNode()) {
                            ++N_Q_NODES;
                            double d = DistanceToNode(e2.GetKey(), e2.GetNodePostfixLen() + 1);
                            N_Q_NODES_0 += d <= 0.0001;
                            if (d <= max_node_dist_) {
                                queue_n_.emplace(d, &e2);
                            }
                        } else {
                            ++N_Q_RESULT;
                            double d = distance_(center_post_, this->post(e2.GetKey()));
                            if (d < max_node_dist_) {
                                queue_v_.emplace(d, &e2);
                                //                                if (queue_v_.size() >
                                //                                num_requested_results_ -
                                //                                num_found_results_) {
                                //                                    // TODO
                                //                                    queue_v_.pop_max();
                                //                                    // TODO do not pop, instead
                                //                                    get max_node_dist from
                                //                                    //
                                //                                    queue_v_[num_requested_results_
                                //                                    - num_found_results_ - 1];
                                //                                }
                                //                                if (queue_v_.size() >=
                                //                                num_requested_results_ -
                                //                                num_found_results_) {
                                //                                    // TODO adjust with 10th value
                                //                                    in queue i.o. last value?
                                //                                    //   -> in case we allow more
                                //                                    than 10...
                                //                                    // TODO !!!!!!!!!!!!!!
                                //                                    max_node_dist_ =
                                //                                    std::min(max_node_dist_,
                                //                                    queue_v_.top_max().first);
                                //                                }
                                if (queue_v_.size() >=
                                    num_requested_results_ - num_found_results_) {
                                    auto pos_max = queue_v_.size() - num_requested_results_ +
                                        num_found_results_;
                                    double d_max = queue_v_[pos_max].first;
                                    max_node_dist_ = std::min(max_node_dist_, d_max);
                                }
                            }
                        }
                    }
                }
                MAX_DEPTH = std::max(MAX_DEPTH, queue_n_.size());
            }
        }
        TOTAL_DEPTH += queue_n_.size();
        this->SetFinished();
        current_distance_ = std::numeric_limits<double>::max();
    }

    double DistanceToNode(const KeyInternal& prefix, std::uint32_t bits_to_ignore) {
        assert(bits_to_ignore < detail::MAX_BIT_WIDTH<SCALAR>);
        SCALAR mask_min = detail::MAX_MASK<SCALAR> << bits_to_ignore;
        SCALAR mask_max = ~mask_min;
        KeyInternal buf;
        // The following calculates the point inside the node that is closest to center_.
        for (dimension_t i = 0; i < DIM; ++i) {
            // if center_[i] is outside the node, return distance to the closest edge,
            // otherwise return center_[i] itself (assume possible distance=0)
            SCALAR min = prefix[i] & mask_min;
            SCALAR max = prefix[i] | mask_max;
            buf[i] = min > center_[i] ? min : (max < center_[i] ? max : center_[i]);
        }
        return distance_(center_post_, this->post(buf));
    }

  private:
    const KeyInternal center_;
    // center after post processing == the external representation
    const KeyExternal center_post_;
    double current_distance_;
    std::
        priority_queue<EntryDistT, std::vector<EntryDistT>, CompareEntryDistByDistance1<EntryDistT>>
            queue_n_;
    //    std::priority_queue<EntryDistT, std::vector<EntryDistT>,
    //    CompareEntryDistByDistance1<EntryDistT>>
    //        queue_v_;
    //    ::phtree::bptree::detail::priority_queue<EntryDistT,
    //    CompareEntryDistByDistance1<EntryDistT>>
    //        queue_n_;
    ::phtree::bptree::detail::priority_queue<EntryDistT, CompareEntryDistByDistance1<EntryDistT>>
        queue_v_;
    size_t num_found_results_;
    size_t num_requested_results_;
    DISTANCE distance_;
    double max_node_dist_ = std::numeric_limits<double>::max();
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_QUERY_KNN_HS1_H

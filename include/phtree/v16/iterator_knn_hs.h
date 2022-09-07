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

#ifndef PHTREE_V16_QUERY_KNN_HS_H
#define PHTREE_V16_QUERY_KNN_HS_H

#include "phtree/common/common.h"
#include "iterator_base.h"
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
using EntryDist = std::pair<double, const Entry<DIM, T, SCALAR>*>;

template <typename ENTRY>
struct CompareEntryDistByDistance {
    bool operator()(const ENTRY& left, const ENTRY& right) const {
        return left.first > right.first;
    };
};
}  // namespace

template <typename T, typename CONVERT, typename DISTANCE, typename FILTER>
class IteratorKnnHS : public IteratorWithFilter<T, CONVERT, FILTER> {
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using KeyExternal = typename CONVERT::KeyExternal;
    using KeyInternal = typename CONVERT::KeyInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = typename IteratorWithFilter<T, CONVERT, FILTER>::EntryT;
    using EntryDistT = EntryDist<DIM, T, SCALAR>;

  public:
    template <typename DIST, typename F>
    explicit IteratorKnnHS(
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
        queue_.emplace(0, &root);
        FindNextElement();
    }

    [[nodiscard]] double distance() const {
        return current_distance_;
    }

    IteratorKnnHS& operator++() noexcept {
        FindNextElement();
        return *this;
    }

    IteratorKnnHS operator++(int) noexcept {
        IteratorKnnHS iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void FindNextElement() {
        while (num_found_results_ < num_requested_results_ && !queue_.empty()) {
            auto& candidate = queue_.top();
            auto* o = candidate.second;
            if (!o->IsNode()) {
                // data entry
                ++num_found_results_;
                this->SetCurrentResult(o);
                current_distance_ = candidate.first;
                // We need to pop() AFTER we processed the value, otherwise the reference is
                // overwritten.
                queue_.pop();
                return;
            } else {
                // inner node
                auto& node = o->GetNode();
                queue_.pop();
                for (auto& entry : node.Entries()) {
                    auto& e2 = entry.second;
                    if (this->ApplyFilter(e2)) {
                        if (e2.IsNode()) {
                            double d = DistanceToNode(e2.GetKey(), e2.GetNodePostfixLen() + 1);
                            queue_.emplace(d, &e2);
                        } else {
                            double d = distance_(center_post_, this->post(e2.GetKey()));
                            queue_.emplace(d, &e2);
                        }
                    }
                }
            }
        }
        this->SetFinished();
        current_distance_ = std::numeric_limits<double>::max();
    }

    double DistanceToNode(const KeyInternal& prefix, std::uint32_t bits_to_ignore) {
        assert(bits_to_ignore < MAX_BIT_WIDTH<SCALAR>);
        SCALAR mask_min = MAX_MASK<SCALAR> << bits_to_ignore;
        SCALAR mask_max = ~mask_min;
        KeyInternal buf;
        // The following calculates the point inside of the node that is closest to center_.
        // If center is inside the node this returns center_, otherwise it finds a point on the
        // node's surface.
        for (dimension_t i = 0; i < DIM; ++i) {
            // if center_[i] is outside the node, return distance to closest edge,
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
    std::priority_queue<EntryDistT, std::vector<EntryDistT>, CompareEntryDistByDistance<EntryDistT>>
        queue_;
    size_t num_found_results_;
    size_t num_requested_results_;
    DISTANCE distance_;
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_QUERY_KNN_HS_H

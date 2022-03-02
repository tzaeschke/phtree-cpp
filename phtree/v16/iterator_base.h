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

#ifndef PHTREE_V16_ITERATOR_BASE_H
#define PHTREE_V16_ITERATOR_BASE_H

#include "../common/common.h"
#include "entry.h"

namespace improbable::phtree::v16 {

template <dimension_t DIM, typename T, typename CONVERT>
class PhTreeV16;

/*
 * Base class for all PH-Tree iterators.
 */
template <typename T, typename CONVERT, typename FILTER = FilterNoOp>
class IteratorBase {
  protected:
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using KeyInternal = typename CONVERT::KeyInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = Entry<DIM, T, SCALAR>;
    friend PhTreeV16<DIM, T, CONVERT>;

  public:
    explicit IteratorBase(const CONVERT& converter)
    : current_result_{nullptr}
    , current_node_{}
    , parent_node_{}
    , is_finished_{false}
    , converter_{converter}
    , filter_{FILTER()} {}

    explicit IteratorBase(const CONVERT& converter, FILTER filter)
    : current_result_{nullptr}
    , current_node_{}
    , parent_node_{}
    , is_finished_{false}
    , converter_{converter}
    , filter_(std::move(filter)) {}

    T& operator*() const {
        assert(current_result_);
        return current_result_->GetValue();
    }

    T* operator->() const {
        assert(current_result_);
        return &current_result_->GetValue();
    }

    template <typename FILTER_2>
    friend bool operator==(
        const IteratorBase<T, CONVERT, FILTER>& left,
        const IteratorBase<T, CONVERT, FILTER_2>& right) {
        // Note: The following compares pointers to Entry objects so it should be
        // a) fast (i.e. not comparing contents of entries)
        // b) return `false` when comparing apparently identical entries from different PH-Trees (as
        // intended)
        return (left.is_finished_ && right.Finished()) ||
            (!left.is_finished_ && !right.Finished() &&
             left.current_result_ == right.GetCurrentResult());
    }

    template <typename FILTER_2>
    friend bool operator!=(
        const IteratorBase<T, CONVERT, FILTER>& left,
        const IteratorBase<T, CONVERT, FILTER_2>& right) {
        return !(left == right);
    }

    auto first() const {
        return converter_.post(current_result_->GetKey());
    }

    T& second() const {
        return current_result_->GetValue();
    }

    [[nodiscard]] bool Finished() const {
        return is_finished_;
    }

    const EntryT* GetCurrentResult() const {
        return current_result_;
    }

  protected:
    void SetFinished() {
        is_finished_ = true;
        current_result_ = nullptr;
    }

    [[nodiscard]] bool ApplyFilter(const EntryT& entry) const {
        return entry.IsNode()
            ? filter_.IsNodeValid(entry.GetKey(), entry.GetNode().GetPostfixLen() + 1)
            : filter_.IsEntryValid(entry.GetKey(), entry.GetValue());
    }

    void SetCurrentResult(const EntryT* current_result) {
        current_result_ = current_result;
    }

    void SetCurrentNodeEntry(const EntryT* current_node) {
        assert(!current_node || current_node->IsNode());
        current_node_ = current_node;
    }

    void SetParentNodeEntry(const EntryT* parent_node) {
        assert(!parent_node || parent_node->IsNode());
        parent_node_ = parent_node;
    }

    auto post(const KeyInternal& point) {
        return converter_.post(point);
    }

  private:
    /*
     * The parent entry contains the parent node. The parent node is the node ABOVE the current node
     * which contains the current entry.
     */
    EntryT* GetCurrentNodeEntry() const {
        return const_cast<EntryT*>(current_node_);
    }

    const EntryT* GetParentNodeEntry() const {
        return parent_node_;
    }

    const EntryT* current_result_;
    const EntryT* current_node_;
    const EntryT* parent_node_;
    bool is_finished_;
    const CONVERT& converter_;
    FILTER filter_;
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ITERATOR_BASE_H

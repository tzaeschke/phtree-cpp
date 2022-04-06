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
template <typename EntryT>
class IteratorBase0 {
    using T = typename EntryT::OrigValueT;

  public:
    explicit IteratorBase0() noexcept : current_result_{nullptr} {}
    explicit IteratorBase0(const EntryT* current_result) noexcept
    : current_result_{current_result} {}

    inline T& operator*() const noexcept {
        assert(current_result_);
        return current_result_->GetValue();
    }

    inline T* operator->() const noexcept {
        assert(current_result_);
        return &current_result_->GetValue();
    }

    inline friend bool operator==(
        const IteratorBase0<EntryT>& left, const IteratorBase0<EntryT>& right) noexcept {
        return left.current_result_ == right.current_result_;
    }

    inline friend bool operator!=(
        const IteratorBase0<EntryT>& left, const IteratorBase0<EntryT>& right) noexcept {
        return left.current_result_ != right.current_result_;
    }

    T& second() const {
        return current_result_->GetValue();
    }

    [[nodiscard]] inline bool IsEnd() const noexcept {
        return current_result_ == nullptr;
    }

    inline const EntryT* GetCurrentResult() const noexcept {
        return current_result_;
    }

  protected:
    void SetFinished() {
        current_result_ = nullptr;
    }

    void SetCurrentResult(const EntryT* current_result) {
        current_result_ = current_result;
    }

  protected:
    const EntryT* current_result_;
};

template <typename EntryT>
using IteratorEnd = IteratorBase0<EntryT>;

template <typename T, typename CONVERT, typename FILTER = FilterNoOp>
class IteratorBase
: public IteratorBase0<Entry<CONVERT::DimInternal, T, typename CONVERT::ScalarInternal>> {
  protected:
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using KeyInternal = typename CONVERT::KeyInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = Entry<DIM, T, SCALAR>;
    friend PhTreeV16<DIM, T, CONVERT>;

  public:
    using FilterT = FILTER;

    explicit IteratorBase(const CONVERT* converter) noexcept
    : IteratorBase0<EntryT>(nullptr)
    , current_node_{}
    , parent_node_{}
    , converter_{converter}
    , filter_{FILTER()} {}

    explicit IteratorBase(const CONVERT* converter, FILTER filter) noexcept
    : IteratorBase0<EntryT>(nullptr)
    , current_node_{}
    , parent_node_{}
    , converter_{converter}
    , filter_(std::forward<FILTER>(filter)) {}

    explicit IteratorBase(
        const EntryT* current_result,
        const EntryT* current_node,
        const EntryT* parent_node,
        const CONVERT* converter) noexcept
    : IteratorBase0<EntryT>(current_result)
    , current_node_{current_node}
    , parent_node_{parent_node}
    , converter_{converter}
    , filter_{FILTER()} {}

    auto first() const {
        return converter_->post(this->current_result_->GetKey());
    }

  protected:
    [[nodiscard]] bool ApplyFilter(const EntryT& entry) const {
        return entry.IsNode() ? filter_.IsNodeValid(entry.GetKey(), entry.GetNodePostfixLen() + 1)
                              : filter_.IsEntryValid(entry.GetKey(), entry.GetValue());
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
        return converter_->post(point);
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

    const EntryT* current_node_;
    const EntryT* parent_node_;
    const CONVERT* converter_;
    FILTER filter_;
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ITERATOR_BASE_H

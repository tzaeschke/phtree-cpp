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

/*
 * Base class for all PH-Tree iterators.
 */
template <typename EntryT>
class IteratorBase {
    using T = typename EntryT::OrigValueT;

  public:
    explicit IteratorBase() noexcept : current_result_{nullptr} {}
    explicit IteratorBase(const EntryT* current_result) noexcept
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
        const IteratorBase<EntryT>& left, const IteratorBase<EntryT>& right) noexcept {
        return left.current_result_ == right.current_result_;
    }

    inline friend bool operator!=(
        const IteratorBase<EntryT>& left, const IteratorBase<EntryT>& right) noexcept {
        return left.current_result_ != right.current_result_;
    }

    T& second() const {
        return current_result_->GetValue();
    }

    [[nodiscard]] inline bool IsEnd() const noexcept {
        return current_result_ == nullptr;
    }

    inline EntryT* GetCurrentResult() const noexcept {
        return const_cast<EntryT*>(current_result_);
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
using IteratorEnd = IteratorBase<EntryT>;

template <typename T, typename CONVERT, typename FILTER = FilterNoOp>
class IteratorWithFilter
: public IteratorBase<Entry<CONVERT::DimInternal, T, typename CONVERT::ScalarInternal>> {
  protected:
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using KeyInternal = typename CONVERT::KeyInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = Entry<DIM, T, SCALAR>;

  public:
    template <typename F>
    explicit IteratorWithFilter(const CONVERT* converter, F&& filter) noexcept
    : IteratorBase<EntryT>(nullptr), converter_{converter}, filter_{std::forward<F>(filter)} {}

    explicit IteratorWithFilter(const EntryT* current_result, const CONVERT* converter) noexcept
    : IteratorBase<EntryT>(current_result), converter_{converter}, filter_{FILTER()} {}

    auto first() const {
        return converter_->post(this->current_result_->GetKey());
    }

    auto& __Filter() {
        return filter_;
    }

  protected:
    [[nodiscard]] bool ApplyFilter(const EntryT& entry) {
        return entry.IsNode() ? filter_.IsNodeValid(entry.GetKey(), entry.GetNodePostfixLen() + 1)
                              : filter_.IsEntryValid(entry.GetKey(), entry.GetValue());
    }

    auto post(const KeyInternal& point) {
        return converter_->post(point);
    }

  private:
    const CONVERT* converter_;
    FILTER filter_;
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ITERATOR_BASE_H

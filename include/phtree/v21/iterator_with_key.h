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

#ifndef PHTREE_V21_ITERATOR_KEY_H
#define PHTREE_V21_ITERATOR_KEY_H

#include "phtree/common/common.h"
#include "iterator_base.h"

namespace improbable::phtree::v20 {

template <dimension_t DIM, typename T, typename CONVERT>
class PhTreeV21;


template <typename T, typename CONVERT>
class IteratorWithKey : public IteratorWithFilter<T, CONVERT> {
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = typename IteratorWithFilter<T, CONVERT>::EntryT;
    using KeyT = typename CONVERT::KeyInternal;
    friend PhTreeV21<DIM, T, CONVERT>; // TODO

  public:
    explicit IteratorWithKey(
        const KeyT& key,
        EntryIteratorC<DIM, EntryT> iter,
        const EntryT* current_result,
        const EntryT* current_node,
        const CONVERT* converter) noexcept
    : IteratorWithFilter<T, CONVERT>(current_result, converter)
    , key_{key}
    , iter_{iter} , current_node_{current_node} {}

    IteratorWithKey& operator++() {
        FindNextElement();
        return *this;
    }

    IteratorWithKey operator++(int) {
        IteratorWithKey iterator(*this);
        ++(*this);
        return iterator;
    }

    EntryT* __GetNodeEntry() const {
        return const_cast<EntryT*>(current_node_);
    }

//    EntryT* __GetParentNodeEntry() const {
//        return const_cast<EntryT*>(parent_node_);
//    }

  private:
    void FindNextElement() noexcept {
        assert(!this->IsEnd());
        auto hc_pos = iter_->first;
        ++iter_;
        while (iter_ != current_node_->GetNode().Entries().end()) {
            if (iter_->first != hc_pos) {
                this->SetFinished();
                return;
            }
            if (iter_->second.GetKey() == key_) {
                this->SetCurrentResult(&iter_->second);
                return;
            }
            ++iter_;
        }
        this->SetFinished();
    }

  private:
    /*
     * The parent entry contains the parent node. The parent node is the node ABOVE the current node
     * which contains the current entry.
     */
    EntryT* GetNodeEntry() const {
        return const_cast<EntryT*>(current_node_);
    }

//    EntryT* GetParentNodeEntry() const {
//        return const_cast<EntryT*>(parent_node_);
//    }

    KeyT key_;
    EntryIteratorC<DIM, EntryT> iter_;
    const EntryT* current_node_;
};

}  // namespace improbable::phtree::v20

#endif  // PHTREE_V21_ITERATOR_KEY_H

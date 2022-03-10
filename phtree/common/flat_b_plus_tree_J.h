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

#ifndef PHTREE_COMMON_B_PLUS_TREE_J_H
#define PHTREE_COMMON_B_PLUS_TREE_J_H

#include "bits.h"
#include <cassert>
#include <tuple>
#include <vector>

/*
 * PLEASE do not include this file directly, it is included via common.h.
 *
 * This file contains the sparse_map implementation, which is used in medium-dimensional nodes in
 * the PH-Tree.
 */
namespace improbable::phtree {

/*
 * The sparse_map is a flat map implementation that uses an array of *at* *most* SIZE=2^DIM.
 * The array contains a list sorted by key.
 *
 * It has O(log n) lookup and O(n) insertion/removal time complexity, space complexity is O(n).
 */

namespace {
template <typename T>
class BptNode;

template <typename T>
struct BptEntry {
    BptNode<T>* node_;
    std::optional<T> value_;
};

template <typename T>
using BptPair = std::pair<size_t, BptEntry<T>>;

using index_t = std::int64_t;

template <typename T>
class b_plus_tree_node {
    using EntryT = BptEntry<T>;

  public:
    explicit b_plus_tree_node(bool is_leaf) : data_{}, is_leaf_{is_leaf} {};

    [[nodiscard]] auto find(size_t key) {
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return it;
        }
        return data_.end();
    }

    [[nodiscard]] auto find(size_t key) const {
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return it;
        }
        return data_.end();
    }

    [[nodiscard]] auto lower_bound(size_t key) {
        return std::lower_bound(
            data_.begin(), data_.end(), key, [](EntryT& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto lower_bound(size_t key) const {
        return std::lower_bound(
            data_.cbegin(), data_.cend(), key, [](const EntryT& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto begin() {
        return data_.begin();
    }

    [[nodiscard]] auto begin() const {
        return cbegin();
    }

    [[nodiscard]] auto cbegin() const {
        return data_.cbegin();
    }

    [[nodiscard]] auto end() {
        return data_.end();
    }

    [[nodiscard]] auto end() const {
        return data_.end();
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return try_emplace_base(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto try_emplace(size_t key, Args&&... args) {
        return try_emplace_base(key, std::forward<Args>(args)...);
    }

    void erase(size_t key) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            data_.erase(it);
        }
    }

    void erase(const typename std::vector<EntryT>::iterator& iterator) {
        data_.erase(iterator);
    }

    [[nodiscard]] size_t size() const {
        return data_.size();
    }

    [[nodiscard]] bool is_leaf() const {
        return is_leaf_;
    }

  private:
    template <typename... Args>
    auto emplace_base(size_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            return std::make_pair(data_.emplace(it, key, std::forward<Args>(args)...), true);
        }
    }

    template <typename... Args>
    auto try_emplace_base(size_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            auto x = data_.emplace(
                it,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
            return std::make_pair(x, true);
        }
    }

    std::vector<EntryT> data_;
    const bool is_leaf_;
};
}  // namespace

constexpr size_t M = 8;

template <typename T>
class b_plus_tree_map {
  public:
    explicit b_plus_tree_map() : root_{true} {};

    [[nodiscard]] auto find(size_t key) {
        auto node = root_;
        while (!node.is_leaf()) {
            auto it = node.find(key);
            if (it == node.end() {
                return end();
            }
            node = it->second;
        }
        auto it = node.lower_bound(key);
        if (it != node.end() && it->first == key) {
            return it;
        }
        return end();
    }

    [[nodiscard]] auto find(size_t key) const {
        auto node = root_;
        while (!node.is_leaf()) {
            auto it = node.find(key);
            if (it == node.end() {
                return end();
            }
            node = it->second;
        }
        auto it = node.lower_bound(key);
        if (it != node.end() && it->first == key) {
            return it;
        }
        return end();
        //        auto it = lower_bound(key);
        //        if (it != data_.end() && it->first == key) {
        //            return it;
        //        }
        //        return data_.end();
    }

    [[nodiscard]] auto lower_bound(size_t key) {
        return std::lower_bound(
            data_.begin(), data_.end(), key, [](PhFlatMapPair<T>& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto lower_bound(size_t key) const {
        return std::lower_bound(
            data_.cbegin(), data_.cend(), key, [](const PhFlatMapPair<T>& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto begin() {
        return data_.begin();
    }

    [[nodiscard]] auto begin() const {
        return cbegin();
    }

    [[nodiscard]] auto cbegin() const {
        return data_.cbegin();
    }

    [[nodiscard]] auto end() {
        return data_.end();
    }

    [[nodiscard]] auto end() const {
        return data_.end();
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return try_emplace_base(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto try_emplace(size_t key, Args&&... args) {
        return try_emplace_base(key, std::forward<Args>(args)...);
    }

    void erase(size_t key) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            data_.erase(it);
        }
    }

    void erase(const typename std::vector<PhFlatMapPair<T>>::iterator& iterator) {
        data_.erase(iterator);
    }

    [[nodiscard]] size_t size() const {
        return data_.size();
    }

  private:
    template <typename... Args>
    auto emplace_base(size_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            return std::make_pair(data_.emplace(it, key, std::forward<Args>(args)...), true);
        }
    }

    template <typename... Args>
    auto try_emplace_base(size_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            auto x = data_.emplace(
                it,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
            return std::make_pair(x, true);
        }
    }

    b_plus_tree_node<T> root_;
};

template <dimension_t DIM>
static std::uint16_t CalcMaxLeafN() {
    // The idea is to have at most one level of inner pages for d<=12
    // The inner pages are all slightly larger the strictly necessary because the fill rate of
    // leaves is < 100%
    switch (DIM) {
    case 1:
        return 2;
    case 2:
        return 4;
    case 3:
        return 8;
    case 4:
        return 16;
    case 5:
        return 16;
    case 6:
        return 16;
    case 7:
        return 16;
    case 8:
        return 16;
    case 9:
        return 32;
    case 10:
        return 32;
    case 11:
        return 32;
    case 12:
        return 64;
    default:
        return 100;
    }
}

template <dimension_t DIM>
static std::uint16_t CalcMaxInnerN() {
    // The idea is to have at most one level of inner pages for d<=12
    // The inner pages are all slightly larger the strictly necessary because the fill rate of
    // leaves is < 100%
    switch (DIM) {
    case 1:
        return 3;
    case 2:
        return 3;
    case 3:
        return 3;
    case 4:
        return 3;
    case 5:
        return 4 + 1;
    case 6:
        return 6 + 1;
    case 7:
        return 10 + 1;
    case 8:
        return 20 + 1;
    case 9:
        return 20 + 1;
    case 10:
        return 35 + 1;
    case 11:
        return 70 + 1;
    case 12:;
        return 70 + 1;
    default:
        return 100;
    }
}

template <dimension_t DIM, typename ENTRY>
class BSTree;

template <dimension_t DIM, typename ENTRY>
class BSTreePage;

template <dimension_t DIM, typename ENTRY>
struct BSTEntry {
    index_t key_;
    ENTRY* value_;
    BSTreePage<DIM, ENTRY>* node_;

    BSTEntry(index_t key, ENTRY* entry) {
        key_ = key;
        kdKey = k;
        value_ = v;
    }

    index_t getKey() {
        return key;
    }

    PhPoint<DIM> getKdKey() {
        return kdKey;
    }

    void* getValue() {
        return value;
    }

    void set(index_t key, PhPoint<DIM> kdKey, void* value) {
        this.key = key;
        this.kdKey = kdKey;
        this.value = value;
    }

    void setValue(void* value) {
        this.value = value;
    }
};

template <dimension_t DIM, typename ENTRY>
class BSTreePage {
    using BSTEntryT = BSTEntry<DIM, ENTRY>;
    using KeyT = PhPoint<DIM>;
    const int INITIAL_PAGE_SIZE = 4;
    const std::uint16_t MAX_LEAF_N = CalcMaxLeafN<DIM>();
    const std::uint16_t MAX_INNER_N = CalcMaxInnerN<DIM>();

    enum REMOVE_OP { REMOVE_RETURN, KEEP_RETURN, KEEP_RETURN_NULL };

    // TODO remove me!
    const int NO_POS = -1000;

    struct UpdateInfo {
        static const int NO_INSERT_REQUIRED = Integer.MAX_VALUE;
        KeyT& newKey;
        int insertRequired = NO_INSERT_REQUIRED;
        UpdateInfo init(KeyT& newKey) {
            this.newKey = newKey;
            return this;
        }
    };

  public:
    BSTreePage(BSTreePage* parent, bool isLeaf, BSTreePage leftPredecessor) {
        init(parent, isLeaf, leftPredecessor);
    }

  private:
    void init(BSTreePage parent, bool isLeaf, BSTreePage leftPredecessor) {
        nextLeaf = nullptr;
        prevLeaf = nullptr;
        this.parent = parent;
        if (isLeaf) {
            nEntries = 0;
            int initialPageSize = MAX_LEAF_N <= 8 ? 2 : INITIAL_PAGE_SIZE;
            keys = tree.bstPool().arrayCreateLong(initialPageSize);
            values = tree.bstPool().arrayCreateEntries(initialPageSize);
            subPages = nullptr;
        } else {
            nEntries = -1;
            keys = tree.bstPool().arrayCreateLong(MAX_INNER_N);
            values = nullptr;
            subPages = tree.bstPool().arrayCreateNodes(MAX_INNER_N + 1);
        }

        this.isLeaf = isLeaf;

        if (isLeaf && leftPredecessor != nullptr) {
            nextLeaf = leftPredecessor.nextLeaf;
            prevLeaf = leftPredecessor;
            leftPredecessor.nextLeaf = this;
            if (nextLeaf != nullptr) {
                nextLeaf.prevLeaf = this;
            }
        } else {
            nextLeaf = nullptr;
            prevLeaf = nullptr;
        }
    }

  public:
    void init(BSTEntryT* e1, BSTEntryT* e2) {
        if (!isLeaf) {
            assert(false && "IllegalStateException");
        }
        if (nEntries > 0) {
            // TODO
            assert(false && "nEntries=" + nEntries);
        }
        values[0] = e1;
        keys[0] = e1.getKey();
        values[1] = e2;
        keys[1] = e2.getKey();
        nEntries = 2;
    }

  public:
    static BSTreePage create(BSTreePage* parent, bool isLeaf, BSTreePage* leftPredecessor) {
        return tree.bstPool().getNode(parent, isLeaf, leftPredecessor);
    }

  public:
    static BSTreePage create(
        BSTreePage* parent, BSTreePage* firstSubpage, BSTreePage* secondSubpage) {
        BSTreePage p = create(parent, false, nullptr);
        p.nEntries++;
        p.subPages[0] = firstSubpage;
        p.nEntries++;
        p.subPages[1] = secondSubpage;
        p.keys[0] = secondSubpage.getMinKey();
        firstSubpage.setParent(p);
        secondSubpage.setParent(p);
        return p;
    }

  private:
    int maxInnerN() {
        return keys.length;
    }

  private:
    int minLeafN(int maxLeafN) {
        return maxLeafN >> 1;
    }

  private:
    int minInnerN(int maxInnerN) {
        return maxInnerN >> 1;
    }

  public:
    BSTreePage* findSubPage(index_t key) {
        // The stored value[i] is the min-values of the according page[i+1}
        int pos = binarySearchInnerNode(key);
        // read page before that value
        return subPages[pos];
    }

  public:
    BSTEntryT* findAndRemove(long key, KeyT& kdKey, Node node, UpdateInfo& ui) {
        // The stored value[i] is the min-values of the according page[i+1}
        int pos = binarySearchInnerNode(key);
        // read page before that value
        BSTreePage* page = subPages[pos];
        BSTEntryT* result;
        if (page->IsLeaf()) {
            result = page->remove(key, kdKey, node, ui);
            checkUnderflowSubpageLeaf(pos, node);
        } else {
            result = page->findAndRemove(key, kdKey, node, ui);
            handleUnderflowSubInner(pos);
        }
        return result;
    }

  public:
    void* getOrCreate(long key, Node ind) {
        // The stored value[i] is the min-values of the according page[i+1}
        int pos = binarySearchInnerNode(key);
        // read page before that value
        BSTreePage page = getPageByPos(pos);
        if (page.IsLeaf()) {
            BSTEntry o = page.getOrCreate(key, this, pos);
            if (o.getKdKey() == nullptr && o.getValue() instanceof BSTreePage) {
                // add page
                BSTreePage newPage = (BSTreePage)o.getValue();
                addSubPage(newPage, newPage.getMinKey(), pos, ind);
                o.setValue(nullptr);
                return o;
            }
            return o;
        }
        return page;
    }

  public:
    BSTEntryT* getValueFromLeaf(long key) {
        int pos = binarySearch(key);
        if (pos >= 0) {
            return values[pos];
        }
        // If the value could is not on this page, it does not exist.
        return nullptr;
    }

  public:
    int binarySearchInnerNode(long key) {
        if (nEntries <= 8) {
            for (int i = 0; i < nEntries; i++) {
                if (key <= keys[i]) {
                    return key == keys[i] ? i + 1 : i;
                }
            }
            return nEntries;  // key not found.
        }
        int low = 0;
        int high = nEntries - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            long midVal = keys[mid];

            if (midVal < key)
                low = mid + 1;
            else if (midVal > key)
                high = mid - 1;
            else {
                return mid + 1;  // key found
            }
        }
        return low;  // key not found.
    }

    /*
     * Binary search.
     *
     * @param key search key
     */
    int binarySearch(long key) {
        if (nEntries <= 8) {
            return linearSearch(key);
        }
        int low = 0;
        int high = nEntries - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            long midVal = keys[mid];

            if (midVal < key)
                low = mid + 1;
            else if (midVal > key)
                high = mid - 1;
            else {
                return mid;  // key found
            }
        }
        return -(low + 1);  // key not found.
    }

  private:
    int linearSearch(long key) {
        for (int i = 0; i < nEntries; i++) {
            if (key <= keys[i]) {
                return key == keys[i] ? i : -(i + 1);
            }
        }
        return -(nEntries + 1);  // key not found.
    }

  private:
    void putUnchecked(int pos, index_t key, BSTEntryT* value) {
        // okay so we add it locally
        shiftArrayForInsertion(pos);
        keys[pos] = key;
        values[pos] = value;
        nEntries++;
        ind.incEntryCount();
    }

  private:
    void shiftArrayForInsertion(int pos) {
        ensureSizePlusOne();
        // Only shift if we do not append
        if (pos < nEntries) {
            System.arraycopy(keys, pos, keys, pos + 1, nEntries - pos);
            System.arraycopy(values, pos, values, pos + 1, nEntries - pos);
        }
    }

  private:
    void ensureSizePlusOne() {
        if (nEntries + 1 > keys.length) {
            int newLen = keys.length * 2 > MAX_LEAF_N ? MAX_LEAF_N : keys.length * 2;
            keys = tree.bstPool().arrayExpand(keys, newLen);
            values = tree.bstPool().arrayExpand(values, newLen);
        }
    }

  private:
    void ensureSize(int newLen) {
        if (newLen > keys.length) {
            keys = tree.bstPool().arrayExpand(keys, newLen);
            values = tree.bstPool().arrayExpand(values, newLen);
        }
    }

  public:
    BSTEntryT* getOrCreate(long key, BSTreePage* parent, int posPageInParent) {
        if (!isLeaf) {
            assert(false && "Tree inconsistency.");
        }

        // in any case, check whether the key(+value) already exists
        int pos = binarySearch(key);
        // key found? -> pos >=0
        if (pos >= 0) {
            return values[pos];
        }
        return create(key, pos, parent, posPageInParent);
    }

  private:
    BSTEntryT* create(index_t key, int pos, BSTreePage* parent, int posPageInParent) {
        BSTEntry value = tree.bstPool().getEntry();
        value.set(key, nullptr, nullptr);

        if (nEntries < MAX_LEAF_N) {
            // okay so we add it locally
            pos = -(pos + 1);
            ensureSizePlusOne();
            if (pos < nEntries) {
                System.arraycopy(keys, pos, keys, pos + 1, nEntries - pos);
                System.arraycopy(values, pos, values, pos + 1, nEntries - pos);
            }
            keys[pos] = key;
            values[pos] = value;
            nEntries++;
            ind.incEntryCount();
            return value;
        }

        // treat page overflow

        // destination page
        BSTreePage destP;
        // created new page?
        bool isNew = false;
        // is previous page?
        bool isPrev = false;

        if (parent == nullptr) {
            destP = ind.bstCreatePage(nullptr, true, this);
            isNew = true;
        } else {
            // use MAX_LEAF_N -1 to avoid pretty much pointless copying (and possible endless
            // loops, see iterator tests)
            BSTreePage next = parent.getNextLeafPage(posPageInParent);
            if (next != nullptr && next.nEntries < MAX_LEAF_N - 1) {
                // merge
                destP = next;
                isPrev = false;
            } else {
                // Merging with prev is not make a big difference, maybe we should remove it...
                BSTreePage prev = parent.getPrevLeafPage(posPageInParent);
                if (prev != nullptr && prev.nEntries < MAX_LEAF_N - 1) {
                    // merge
                    destP = prev;
                    isPrev = true;
                } else {
                    destP = ind.bstCreatePage(parent, true, this);
                    isNew = true;
                }
            }
        }

        // Ensure all nodes have full capacity
        ensureSize(MAX_LEAF_N);
        destP.ensureSize(MAX_LEAF_N);

        // We move 50% of data. For bulkloading, we could keep 95% or so in old page. 100%? But
        // there is no bulk loading.
        int nEntriesToKeep = (nEntries + destP.nEntries) >> 1;
        int nEntriesToCopy = nEntries - nEntriesToKeep;
        if (isNew) {
            // works only if new page follows current page
            System.arraycopy(keys, nEntriesToKeep, destP.keys, 0, nEntriesToCopy);
            System.arraycopy(values, nEntriesToKeep, destP.values, 0, nEntriesToCopy);
        } else if (isPrev) {
            // copy element to previous page
            System.arraycopy(keys, 0, destP.keys, destP.nEntries, nEntriesToCopy);
            System.arraycopy(values, 0, destP.values, destP.nEntries, nEntriesToCopy);
            // move element forward to beginning of page
            System.arraycopy(keys, nEntriesToCopy, keys, 0, nEntries - nEntriesToCopy);
            System.arraycopy(values, nEntriesToCopy, values, 0, nEntries - nEntriesToCopy);
        } else {
            // make space on next page
            System.arraycopy(destP.keys, 0, destP.keys, nEntriesToCopy, destP.nEntries);
            System.arraycopy(destP.values, 0, destP.values, nEntriesToCopy, destP.nEntries);
            // insert element in next page
            System.arraycopy(keys, nEntriesToKeep, destP.keys, 0, nEntriesToCopy);
            System.arraycopy(values, nEntriesToKeep, destP.values, 0, nEntriesToCopy);
        }
        pos = -(pos + 1);
        int oldNEntriesP = destP.nEntries;
        nEntries = nEntriesToKeep;
        destP.nEntries = nEntriesToCopy + destP.nEntries;
        // New page and min key
        if (isNew || !isPrev) {
            if (destP.keys[0] > key) {
                putUnchecked(pos, key, value, ind);
            } else {
                destP.putUnchecked(pos - nEntriesToKeep, key, value, ind);
            }
        } else {
            if (keys[0] > key) {
                destP.putUnchecked(pos + oldNEntriesP, key, value, ind);
            } else {
                putUnchecked(pos - nEntriesToCopy, key, value, ind);
            }
        }
        if (isNew) {
            // own key remains unchanged
            // Hack: we return the new page as value of BSEntry
            value.setValue(destP);
            return value;
        } else {
            // change own key in parent?
            if (isPrev) {
                // posInParent has not changed!
                parent.updateKey(getMinKey(), posPageInParent - 1);
            } else {
                // change key of 'next' page
                parent.updateKey(destP.getMinKey(), posPageInParent - 1 + 1);
            }
            return value;
        }
    }

  private:
    void updateKey(long key, int keyPos) {
        if (keyPos < 0 || keys[keyPos] == key) {
            // nothing changes
            return;
        }

        keys[keyPos] = key;
    }

  private:
    void addSubPage(BSTreePage newP, long minKey, int keyPos, Node ind) {
        if (isLeaf) {
            assert(false && "Tree inconsistency");
        }

        if (nEntries < maxInnerN()) {
            // add page here

            if (keyPos == NO_POS) {
                // For now, we assume a unique index.
                keyPos = binarySearchInnerNode(minKey);
                // If the key has a perfect match then something went wrong. This should
                // never happen so we don't need to check whether (i < 0).
            }

            if (keyPos > 0) {
                System.arraycopy(keys, keyPos, keys, keyPos + 1, nEntries - keyPos);
                System.arraycopy(subPages, keyPos + 1, subPages, keyPos + 2, nEntries - keyPos);
                keys[keyPos] = minKey;
                subPages[keyPos + 1] = newP;
                newP.setParent(this);
                nEntries++;
            } else {
                // decide whether before or after first page (both will end up before the current
                // first key).
                int ii;
                if (nEntries < 0) {
                    // can happen for empty root page
                    ii = 0;
                } else {
                    System.arraycopy(keys, 0, keys, 1, nEntries);
                    long oldKey = subPages[0].getMinKey();
                    if (minKey > oldKey) {
                        ii = 1;
                        keys[0] = minKey;
                    } else {
                        ii = 0;
                        keys[0] = oldKey;
                    }
                    System.arraycopy(subPages, ii, subPages, ii + 1, nEntries - ii + 1);
                }
                subPages[ii] = newP;
                newP.setParent(this);
                nEntries++;
            }
        } else {
            // treat page overflow
            BSTreePage* newInner = ind.bstCreatePage(parent, false, nullptr);

            // TODO use optimized fill ratio for unique values, just like for leaves?.
            int minInnerN = minInnerN(keys.length);
            System.arraycopy(keys, minInnerN + 1, newInner.keys, 0, nEntries - minInnerN - 1);
            System.arraycopy(subPages, minInnerN + 1, newInner.subPages, 0, nEntries - minInnerN);
            newInner->nEntries = (short)(nEntries - minInnerN - 1);
            newInner->assignThisAsParentToLeaves();

            if (parent == nullptr) {
                // create a parent
                BSTreePage newRoot = ind.bstCreatePage(nullptr, false, nullptr);
                newRoot.subPages[0] = this;
                newRoot.nEntries = 0;  // 0: indicates one leaf / zero keys
                setParent(newRoot);
                ind.bstUpdateRoot(newRoot);
            }

            parent->addSubPage(newInner, keys[minInnerN], NO_POS, ind);

            nEntries = (short)(minInnerN);
            // finally add the leaf to the according page
            BSTreePage* newHome;
            long newInnerMinKey = newInner.getMinKey();
            if (minKey < newInnerMinKey) {
                newHome = this;
            } else {
                newHome = newInner;
            }
            newHome->addSubPage(newP, minKey, NO_POS, ind);
        }
    }

  private:
    long getMinKey() {
        if (isLeaf) {
            return keys[0];
        }
        return getPageByPos(0).getMinKey();
    }

  public:
    void print(String indent) {
        if (isLeaf) {
            System.out.println(
                indent + "Leaf page: nK=" + nEntries + " keys=" + Arrays.toString(keys));
            System.out.println(indent + "                         " + Arrays.toString(values));
        } else {
            System.out.println(
                indent + "Inner page: nK=" + nEntries + " keys=" + Arrays.toString(keys));
            //			System.out.println(indent + "                " + nEntries + " leaf=" +
            //					Arrays.toString(leaves));
            System.out.print(indent + "[");
            for (int i = 0; i <= nEntries; i++) {
                if (subPages[i] != nullptr) {
                    System.out.print(indent + "i=" + i + ": ");
                    subPages[i].print(indent + "  ");
                }
            }
            System.out.println(']');
        }
    }

  public:
    void toStringTree(StringBuilderLn sb, String indent) {
        if (isLeaf) {
            sb.appendLn(indent + "Leaf page: nK=" + nEntries + " keys=" + Arrays.toString(keys));
            sb.appendLn(indent + "                         " + Arrays.toString(values));
        } else {
            sb.appendLn(indent + "Inner page: nK=" + nEntries + " keys=" + Arrays.toString(keys));
            //			System.out.println(indent + "                " + nEntries + " leaf=" +
            //					Arrays.toString(leaves));
            sb.append(indent + "[");
            for (int i = 0; i <= nEntries; i++) {
                if (subPages[i] != nullptr) {
                    sb.append(indent + "i=" + i + ": ");
                    subPages[i].toStringTree(sb, indent + "  ");
                }
            }
            sb.appendLn("]");
        }
    }

  public:
    void printLocal() {
        System.out.println("PrintLocal() for " + this);
        if (isLeaf) {
            System.out.println("Leaf page: nK=" + nEntries + " oids=" + Arrays.toString(keys));
            System.out.println("                         " + Arrays.toString(values));
        } else {
            System.out.println("Inner page: nK=" + nEntries + " oids=" + Arrays.toString(keys));
            System.out.println("                      " + Arrays.toString(subPages));
        }
    }

  public:
    short getNKeys() {
        return nEntries;
    }

  public:
    BSTEntryT* remove(index_t key, KeyT& kdKey, Node node, UpdateInfo& ui) {
        int i = binarySearch(key);
        if (i < 0) {
            // key not found
            return nullptr;
        }

        // first remove the element
        BSTEntryT* prevValue = values[i];
        REMOVE_OP op = node.bstInternalRemoveCallback(prevValue, kdKey, ui);
        switch (op) {
        case REMOVE_RETURN:
            System.arraycopy(keys, i + 1, keys, i, nEntries - i - 1);
            System.arraycopy(values, i + 1, values, i, nEntries - i - 1);
            nEntries--;
            node.decEntryCount();
            return prevValue;
        case KEEP_RETURN:
            return prevValue;
        case KEEP_RETURN_NULL:
            return nullptr;
        default:
            assert(false && "IllegalArgumentException");
        }
    }

  public void* computeLeaf(long key, KeyT& kdKey, int posInParent, Node node,
                                    bool doIfAbsent, BiFunction<long[], ? super T, ? extends T> mappingFunction) {
        int pos = binarySearch(key);
        if (pos < 0) {
            // key not found
            if (doIfAbsent) {
                T newValue = mappingFunction.apply(kdKey, nullptr);
                if (newValue != nullptr) {
                    BSTEntry e = addForCompute(key, pos, posInParent, node);
                    e.set(key, kdKey, newValue);
                    return newValue;
                }
            }
            return null;
        }

        BSTEntryT* currentEntry = values[pos];
        void* currentValue = currentEntry.getValue();
        if (currentValue instanceof Node) {
            if (((Node)currentValue).getInfixLen() == 0) {
                // Shortcut that avoid MCB calculation: No infix conflict, just traverse the subnode
                // (=currentValue)
                return currentValue;
            }
        }

        long[] localKdKey = currentEntry.getKdKey();
        int maxConflictingBits = Node.calcConflictingBits(kdKey, localKdKey);
        if (maxConflictingBits == 0) {
            if (currentValue instanceof Node) {
                // return entry with subnode
                return currentValue;
            }
            T newValue =
                mappingFunction.apply(kdKey, PhTreeHelper.unmaskNull(currentEntry.getValue()));
            if (newValue == nullptr) {
                // remove
                removeForCompute(key, pos, posInParent, node);
                tree.bstPool().offerEntry(currentEntry);
                return nullptr;
            } else {
                // replace (cannot be null)
                currentEntry.setValue(newValue);
            }
            return newValue;
        }

        if (currentValue instanceof Node) {
            Node subNode = (Node)currentValue;
            if (subNode.getPostLen() + 1 >= maxConflictingBits) {
                return subNode;
            }
        }

        // Key found, but entry does not match
        if (doIfAbsent) {
            // We have two entries in the same location (local hcPos).
            // If the kdKey differs, we have to split, insert a newSubNode and return null.
            T newValue = mappingFunction.apply(kdKey, nullptr);
            if (newValue != nullptr) {
                insertSplit(currentEntry, kdKey, newValue, tree, maxConflictingBits, node);
                return newValue;
            }
            return nullptr;
        }
        // Return 'null' when ignoring absent values
        return nullptr;
    }

  private:
    BSTEntryT* addForCompute(long key, int pos, int posPageInParent, Node node) {
        BSTEntryT* o = create(key, pos, parent, posPageInParent, node);
        Node.incEntryCountTree(tree);
        if (o.getKdKey() == nullptr && o.getValue() instanceof BSTreePage) {
            // add page
            BSTreePage newPage = (BSTreePage)o.getValue();
            if (parent != nullptr) {
                parent.addSubPage(newPage, newPage.getMinKey(), posPageInParent, node);
            } else {
                node.bstSetRoot(create(node, nullptr, this, newPage, tree));
            }
            o.setValue(nullptr);
            return o;
        }
        return o;
    }

  private:
    void removeForCompute(long key, int pos, int posPageInParent, Node node) {
        int i = pos;
        System.arraycopy(keys, i + 1, keys, i, nEntries - i - 1);
        System.arraycopy(values, i + 1, values, i, nEntries - i - 1);
        nEntries--;
        node.decEntryCountGlobal(tree);
        if (parent == nullptr) {
            return;
        }
        BSTreePage parentPage = parent;
        parentPage.checkUnderflowSubpageLeaf(posPageInParent, node);

        parentPage = parentPage.parent;
        while (parentPage != nullptr) {
            pos = parentPage.binarySearchInnerNode(key);
            parentPage.handleUnderflowSubInner(pos);
            parentPage = parentPage.parent;
        }
    }

  private:
    void insertSplit(
        BSTEntry currentEntry,
        const KeyT& newKey,
        void* newValue,
        int maxConflictingBits,
        Node node) {
        KeyT& localKdKey = currentEntry.getKdKey();
        Node newNode = node.createNode(
            newKey, newValue, localKdKey, currentEntry.getValue(), maxConflictingBits, tree);
        // replace local entry with new subnode
        currentEntry.set(currentEntry.getKey(), tree.longPool().arrayClone(localKdKey), newNode);
        Node.incEntryCountTree(tree);
    }

  private:
    void checkUnderflowSubpageLeaf(int pos) {
        BSTreePage subPage = getPageByPos(pos);
        if (subPage.nEntries == 0) {
            removePage(pos);
        } else if (subPage.nEntries < minLeafN(MAX_LEAF_N) && (subPage.nEntries % 8 == 0)) {
            // The second term prevents frequent reading of previous and following pages.
            // TODO Should we instead check for nEntries==MAx>>1 then == (MAX>>2) then <= (MAX>>3)?

            // now attempt merging this page
            BSTreePage prevPage = getPrevLeafPage(pos);
            if (prevPage != nullptr) {
                // We merge only if they all fit on a single page. This means we may read
                // the previous page unnecessarily, but we avoid writing it as long as
                // possible. TODO find a balance, and do no read prev page in all cases
                if (subPage.nEntries + prevPage.nEntries < MAX_LEAF_N) {
                    // TODO for now this work only for leaves with the same root. We
                    // would need to update the min values in the inner nodes.
                    System.arraycopy(
                        subPage.keys, 0, prevPage.keys, prevPage.nEntries, subPage.nEntries);
                    System.arraycopy(
                        subPage.values, 0, prevPage.values, prevPage.nEntries, subPage.nEntries);
                    prevPage.nEntries += subPage.nEntries;
                    removePage(pos);
                }
            }
        }
    }

  private:
    void removePage(int posToRemove) {
        BSTreePage indexPage = getPageByPos(posToRemove);

        // remove sub page page from FSM.
        tree.bstPool().reportFreeNode(indexPage);

        if (nEntries > 0) {  // otherwise we just delete this page
            // remove entry
            arraysRemoveInnerEntry(posToRemove);
            nEntries--;
        }
    }

  private:
    void handleUnderflowSubInner(int pos) {
        BSTreePage* sub = getPageByPos(pos);
        if (sub.nEntries < maxInnerN() >> 1) {
            if (sub.nEntries >= 0) {
                BSTreePage prev = getPrevInnerPage(pos);
                if (prev != nullptr && !prev.isLeaf) {
                    // this is only good for merging inside the same parent.
                    if ((sub.nEntries % 2 == 0) && (prev.nEntries + sub.nEntries < maxInnerN())) {
                        System.arraycopy(sub.keys, 0, prev.keys, prev.nEntries + 1, sub.nEntries);
                        System.arraycopy(
                            sub.subPages, 0, prev.subPages, prev.nEntries + 1, sub.nEntries + 1);
                        // find key for the first appended page -> go up or go down????? Up!
                        prev.keys[prev.nEntries] = keys[pos - 1];
                        prev.nEntries += sub.nEntries + 1;  // for the additional key
                        prev.assignThisAsParentToLeaves();
                        removePage(pos);
                    }
                    return;
                }

                if (sub.nEntries == 0) {
                    // only one element left, no merging occurred -> move sub-page up to parent
                    BSTreePage child = sub.getPageByPos(0);
                    replaceChildPage(child, pos);
                    tree.bstPool().reportFreeNode(sub);
                }
            } else {
                // nEntries == 0
                if (sub.parent != nullptr) {
                    return;
                }
                // else : No root and this is a leaf page... -> we do nothing.
                sub.subPages[0] = nullptr;
                sub.nEntries--;  // down to -1 which indicates an empty root page
            }
        }
    }

  private:
    void arraysRemoveKey(int pos) {
        System.arraycopy(keys, pos + 1, keys, pos, nEntries - pos - 1);
    }

  private:
    void arraysRemoveChild(int pos) {
        System.arraycopy(subPages, pos + 1, subPages, pos, nEntries - pos);
        subPages[nEntries] = null;
    }

    /**
     *
     * @param posEntry The pos in the subPage-array. The according keyPos may be -1.
     */
  private:
    void arraysRemoveInnerEntry(int posEntry) {
        if (posEntry > 0) {
            arraysRemoveKey(posEntry - 1);
        } else {
            arraysRemoveKey(0);
        }
        arraysRemoveChild(posEntry);
    }

    /**
     * Replacing sub-pages occurs when the sub-page shrinks down to a single sub-sub-page, in which
     * case we pull up the sub-sub-page to the local page, replacing the sub-page.
     */
  private:
    void replaceChildPage(BSTreePage subChild, int pos) {
        subPages[pos] = subChild;
        if (pos > 0) {
            keys[pos - 1] = subChild.getMinKey();
        }
        subChild.setParent(this);
    }

  public:
    void setParent(BSTreePage parent) {
        this.parent = parent;
    }

    KeyT& getKeys() {
        return keys;
    }

    BSTEntryT[] getValues() {
        return values;
    }

  private:
    void setNEntries(int n) {
        nEntries = (short)n;
    }

    /**
     * Returns only INNER pages.
     * TODO for now this ignores leafPages on a previous inner node. It returns only leaf pages
     * from the current node.
     * @param currentSubPos pos in subnode
     * @return The previous leaf page or null, if the given page is the first page.
     */
  private:
    BSTreePage* getPrevInnerPage(int currentSubPos) {
        if (currentSubPos > 0) {
            BSTreePage page = getPageByPos(currentSubPos - 1);
            if (page.isLeaf) {
                return nullptr;
            }
            return page;
        }
        // TODO we really should return the last leaf page of the previous inner page.
        // But if they get merged, then we have to shift minimal values, which is
        // complicated. For now we ignore this case, hoping that it doesn't matter much.
        return nullptr;
    }

    /**
     * Returns only LEAF pages.
     * @param currentSubPos pos in subnode
     * @return The previous leaf page or null, if the given page is the first page.
     */
  private:
    BSTreePage* getPrevLeafPage(int currentSubPos) {
        if (currentSubPos > 0) {
            BSTreePage page = getPageByPos(currentSubPos - 1);
            return page.getLastLeafPage();
        }
        return nullptr;
    }

    /**
     * Returns only LEAF pages.
     * @param currentSubPos pos in subnode
     * @return The previous next page or null, if the given page is the first page.
     */
  private:
    BSTreePage* getNextLeafPage(int currentSubPos) {
        // do not add +1, because we compare pages to keys
        if (currentSubPos < getNKeys()) {
            return getPageByPos(currentSubPos + 1).getFirstLeafPage();
        }
        return nullptr;
    }

    /**
     *
     * @return The first leaf page of this branch.
     */
  private:
    BSTreePage* getFirstLeafPage() {
        if (isLeaf) {
            return this;
        }
        return getPageByPos(0).getFirstLeafPage();
    }

    /**
     *
     * @return The last leaf page of this branch.
     */
  private:
    BSTreePage* getLastLeafPage() {
        if (isLeaf) {
            return this;
        }
        return getPageByPos(getNKeys()).getLastLeafPage();
    }

    /**
     * Returns the page at the specified position.
     */
    BSTreePage getPageByPos(int pos) {
        return subPages[pos];
    }

  private:
    void assignThisAsParentToLeaves() {
        for (int i = 0; i <= getNKeys(); i++) {
            // leaves may be null if they are not loaded!
            if (subPages[i] != nullptr) {
                subPages[i].setParent(this);
            }
        }
    }

  public:
    void clear() {
        if (!isLeaf) {
            for (int i = 0; i < getNKeys() + 1; i++) {
                BSTreePage p = getPageByPos(i);
                p.clear();
                // 0-IDs are automatically ignored.
                tree.bstPool().reportFreeNode(p);
            }
        }
        if (subPages != nullptr) {
            for (int i = 0; i < subPages.length; i++) {
                subPages[i] = nullptr;
            }
        }
        setNEntries(-1);
    }

  public:
    bool IsLeaf() {
        return isLeaf;
    }

  public:
    void getStats(BSTStats stats) {
        if (IsLeaf()) {
            stats.nNodesLeaf++;
            stats.nEntriesLeaf += nEntries;
            stats.capacityLeaf += keys.length;
            if (nEntries < 1 && parent != nullptr && parent->parent != nullptr) {
                assert(false && "IllegalStateException");
            }
        } else {
            stats.nNodesInner++;
            stats.nEntriesInner += nEntries + 1;
            stats.capacityInner += keys.length + 1;
            if (nEntries < 1 && parent != nullptr) {
                assert(false && "IllegalStateException");
            }
            for (int i = 0; i < getNKeys() + 1; i++) {
                getPageByPos(i).getStats(stats);
            }
        }
    }

  public
    BSTEntry getFirstValue() {
        return values[0];
    }

  public
    BSTreePage getFirstSubPage() {
        return subPages[0];
    }

  public
    BSTreePage[] getSubPages() {
        return subPages;
    }

    void nullify() {
        keys = nullptr;
        values = nullptr;
        subPages = nullptr;
        nextLeaf = nullptr;
        prevLeaf = nullptr;
        parent = nullptr;
        nEntries = 0;
    }

    BSTreePage getNextLeaf() {
        return nextLeaf;
    }

    void updateNeighborsRemove() {
        if (prevLeaf != nullptr) {
            prevLeaf.nextLeaf = nextLeaf;
        }
        if (nextLeaf != nullptr) {
            nextLeaf.prevLeaf = prevLeaf;
        }
    }

  private:
    BSTreePage* parent;
    std::vector<index_t> keys;
    std::vector<BSTEntryT*> values;
    /** number of keys. There are nEntries+1 subPages in any leaf page. */
    size_t nEntries;

    bool isLeaf;
    std::vector<BSTreePage*> subPages;
    BSTreePage* prevLeaf;
    BSTreePage* nextLeaf;
};

template <dimension_t DIM, typename ENTRY>
class BSTree {
    using PageT = BSTreePage<DIM, ENTRY>;

  public:
    BSTree() : root_{new BSTreePage<DIM, ENTRY>()}, size_{0} {}

  private:
    PageT* root_;
    size_t size_;
};

};  // namespace improbable::phtree

#endif  // PHTREE_COMMON_B_PLUS_TREE_J_H

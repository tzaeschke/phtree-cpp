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

#ifndef PHTREE_COMMON_BASE_TYPES_H
#define PHTREE_COMMON_BASE_TYPES_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <sstream>
#include <vector>

/*
 * PLEASE do not include this file directly, it is included via common.h.
 *
 * This file contains specifications for various types used in the PH-Tree, including
 * PhPoint, PhPointD and PhPointBox.
 */
namespace improbable::phtree {

// ************************************************************************
// SFINAE helper
// ************************************************************************
template <typename T>
constexpr auto hc_pos_t_helper(int)
    -> decltype(std::integral_constant<decltype(T::size()), T::size()>{}, typename std::conditional_t<T{}.size() < 32, std::uint32_t, std::uint64_t>{});
template <typename>
constexpr auto hc_pos_t_helper(long) -> std::uint64_t;
template <typename T>
// This returns 32bit for point types with `constexpr size()` and DIM < 32, otherwise 64bit.
using hc_pos_point_t = decltype(hc_pos_t_helper<T>(0));

// ************************************************************************
// Constants and base types
// ************************************************************************

using scalar_64_t = int64_t;
using scalar_32_t = int32_t;
using scalar_16_t = int16_t;

// Bits in a coordinate (usually a double or long has 64 bits, so uint_8 suffices).
// However, uint32_t turned out to be faster, probably due to fewer cycles required for 32bit
// instructions (8bit/16bit tend to require more cycles, see CPU tables available on the web).
using bit_width_t = uint32_t;
// Number of bit for 'scalar_64_t' or 'scalar_32_t'. Note that 'digits' does _not_ include sign bit,
// so e.g. int64_t has 63 `digits`, however we need all bits, i.e. 64.
template <typename SCALAR>
static constexpr bit_width_t MAX_BIT_WIDTH =
    std::numeric_limits<SCALAR>::digits + std::numeric_limits<SCALAR>::is_signed;
// Bit mask
template <typename SCALAR>
using bit_mask_t = typename std::make_unsigned<SCALAR>::type;
template <typename SCALAR>
static constexpr bit_mask_t<SCALAR> MAX_MASK = std::numeric_limits<bit_mask_t<SCALAR>>::max();
using dimension_t = size_t;  // Number of dimensions
// We have two types that represent hypercube addresses (HC position).
// The hc_pos_dim_t uses a template parameter to determine how many bits are needed, this is either
// 32bit or 64bit. This parameter is used where HC positions are stored because benchmarks show a
// difference in performance when this is used.
// The hc_pos_64_t type is always set to 64. It is used where computations play a role that appear
// to prefer being in always 64bit, mainly in CalcPosInArray() and in Node.
template <dimension_t DIM>
using hc_pos_dim_t = std::conditional_t<(DIM < 32), uint32_t, uint64_t>;
using hc_pos_64_t = uint64_t;
// template <typename POINT>
// using hc_pos_point_t = std::conditional_t<std::is_same_v<isConstExprSize<POINT>, std::true_type>,
// hc_pos_dim_t<POINT{}.size()>, uint64_t> ;

// ************************************************************************
// Basic structs and classes
// ************************************************************************

template <dimension_t DIM, typename SCALAR = scalar_64_t>
struct PhPointLD : public std::array<SCALAR, DIM> {};
template <dimension_t DIM, typename SCALAR = scalar_64_t>
// struct PhPointHD : public std::vector<SCALAR> {};
struct PhPointHD {
    using value_type = SCALAR;

  public:
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator = value_type*;
    using const_iterator = const value_type*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    // typedef std::reverse_iterator<iterator>	      reverse_iterator;
    // typedef std::reverse_iterator<const_iterator>   const_reverse_iterator;

    PhPointHD() noexcept {
        array_ = new std::array<SCALAR, DIM>;
    };

    PhPointHD(const PhPointHD& orig) noexcept {
        array_ = new std::array<SCALAR, DIM>(*orig.array_);
    };

    //    PhPointHD(const std::array<SCALAR, DIM>& array) noexcept {
    //        array_ = new std::array<SCALAR, DIM>(array);
    //    };

    PhPointHD(PhPointHD&& orig) noexcept {
        array_ = orig.array_;
        orig.array_ = nullptr;  // TODO remove
    };

    // PhBox(const std::array<SCALAR, DIM>& min, const std::array<SCALAR, DIM>& max)
    //: min_{min}, max_{max} {} // TODO put back in
    // TODO decide:
    //  - either use inheritance of PhPointLD:
    //      -> Can be used as std:array externally
    //  - use delegation
    //      -> PhBox and other constructors function continue to accept std::array as argument....

//    explicit PhPointHD(const std::array<SCALAR, DIM>& array) noexcept {
//        array_ = new std::array<SCALAR, DIM>(array);
//    }

//    template<typename ...E>
//    PhPointHD(E&&...e) : array_{new std::array<SCALAR, DIM>{{std::forward<E>(e)...}}} {}

    // TODO this ALMOST works.....
//    template<typename E1, typename ...E>
//    PhPointHD(E1 e1, E...e) : array_{new std::array<SCALAR, DIM>{{e1, e...}}} {}

    //    template<typename ...E>
//    PhPointHD(E...e) : array_{new std::array<SCALAR, DIM>{{e...}}} {}


    //template<typename E1, typename E2, typename ...E>
    PhPointHD(SCALAR e) : array_{new std::array<SCALAR, DIM>{{e}}} {
        // TODO enable_if

    }

    template<typename E1, typename E2, typename ...E>
    PhPointHD(E1 e1, E2 e2, E...e) : array_{new std::array<SCALAR, DIM>{{e1, e2, e...}}} {
       // TODO enable_if
        static_assert(std::is_same_v<E1, E2>);
    }


//    template <typename... Args>
//    explicit PhPointHD(Args&&... t) : array_{new std::array<SCALAR, DIM>{{std::forward<Args>(t)...}}} {}

//    PhPointHD(std::initializer_list<SCALAR> inputs)
//    :array_{new std::array<SCALAR, DIM>{inputs}} {}

    PhPointHD& operator=(const PhPointHD<DIM, SCALAR>& other) noexcept {
        if (this == &other) {
            return *this;
        }
        *array_ = *other.array_;
        return *this;
    }

    PhPointHD& operator=(PhPointHD<DIM, SCALAR>&& other) noexcept {
        array_ = other.array_;
        other.array_ = nullptr;  // TODO remove
        return *this;
    }

    ~PhPointHD() {
        delete array_;
    }

    constexpr SCALAR& operator[](size_t index) noexcept {
        return (*array_)[index];
    }

    const constexpr SCALAR& operator[](size_t index) const noexcept {
        return (*array_)[index];
    }

    auto operator==(const PhPointHD<DIM, SCALAR>& other) const noexcept {
        return *array_ == *other.array_;
    }

    auto operator!=(const PhPointHD<DIM, SCALAR>& other) const noexcept {
        return *array_ != *other.array_;
    }

    constexpr bool operator<(const PhPointHD& other) const noexcept {
        return *array_ < *other.array_;
    }

    constexpr iterator begin() noexcept {
        return iterator(array_->begin());
    }

    constexpr const_iterator begin() const noexcept {
        return const_iterator(array_->begin());
    }

    constexpr iterator end() noexcept {
        return iterator(array_->end());
    }

    constexpr const_iterator end() const noexcept {
        return const_iterator(array_->end());
    }

    [[nodiscard]] constexpr size_t size() const noexcept {
        return DIM;
    }

  private:
    std::array<SCALAR, DIM>* array_;
};
// template<typename SCALAR, typename... _Up>
// PhPointLD(SCALAR, _Up...)
//     -> PhPointLD<std::enable_if_t<(std::is_same_v<SCALAR, _Up> && ...), SCALAR>,
//                  1 + sizeof...(_Up)>;

// The SCALAR type needs to be a signet integer, i.e. int32_t or int64_t.
template <dimension_t DIM, typename SCALAR = scalar_64_t>
// using PhPoint = std::conditional_t < DIM<300, PhPointLD<DIM, SCALAR>, PhPointHD<DIM, SCALAR>>;
//   using PhPoint = std::array<SCALAR, DIM>;
using PhPoint = PhPointHD<DIM, SCALAR>;

template <dimension_t DIM>
using PhPointD = PhPoint<DIM, double>;

template <dimension_t DIM>
using PhPointF = PhPoint<DIM, float>;

template <dimension_t DIM, typename SCALAR = scalar_64_t>
class PhBox {
    using Point = PhPoint<DIM, SCALAR>;

  public:
    explicit PhBox() = default;

    PhBox(const PhBox<DIM, SCALAR>& orig) = default;

    // PhBox(const std::array<SCALAR, DIM>& min, const std::array<SCALAR, DIM>& max)
    //: min_{min}, max_{max} {} // TODO put back in
    // TODO decide:
    //  - either use inheritance of PhPointLD:
    //      -> Can be used as std:array externally
    //  - use delegation
    //      -> PhBox and other constructors function continue to accept std::array as argument....

//    PhBox(const PhPointLD<DIM, SCALAR>& min, const PhPointLD<DIM, SCALAR>& max)
//    : min_{min}, max_{max} {}
//
//    // TODO avoid copy??
//    PhBox(const PhPointHD<DIM, SCALAR>& min, const PhPointHD<DIM, SCALAR>& max)
//    : min_{min}, max_{max} {}

    PhBox(const PhPoint<DIM, SCALAR>& min, const PhPoint<DIM, SCALAR>& max)
    : min_{min}, max_{max} {}

    [[nodiscard]] const Point& min() const {
        return min_;
    }

    [[nodiscard]] const Point& max() const {
        return max_;
    }

    [[nodiscard]] Point& min() {
        return min_;
    }

    [[nodiscard]] Point& max() {
        return max_;
    }

    void min(const std::array<SCALAR, DIM>& new_min) {
        min_ = new_min;
    }

    void max(const std::array<SCALAR, DIM>& new_max) {
        max_ = new_max;
    }

    auto operator==(const PhBox<DIM, SCALAR>& other) const -> bool {
        return min_ == other.min_ && max_ == other.max_;
    }

    auto operator!=(const PhBox<DIM, SCALAR>& other) const -> bool {
        return !(*this == other);
    }

  private:
    Point min_;
    Point max_;
};

template <dimension_t DIM>
using PhBoxD = PhBox<DIM, double>;

template <dimension_t DIM>
using PhBoxF = PhBox<DIM, float>;

template <improbable::phtree::dimension_t DIM, typename SCALAR>
std::ostream& operator<<(std::ostream& os, const PhPointLD<DIM, SCALAR>& data) {
    assert(data.size() >= 1);
    os << "[";
    for (dimension_t i = 0; i < data.size() - 1; ++i) {
        os << data[i] << ",";
    }
    os << data[data.size() - 1] << "]";
    return os;
}

template <improbable::phtree::dimension_t DIM, typename SCALAR>
std::ostream& operator<<(std::ostream& os, const PhPointHD<DIM, SCALAR>& data) {
    assert(data.size() >= 1);
    os << "[";
    for (dimension_t i = 0; i < data.size() - 1; ++i) {
        os << data[i] << ",";
    }
    os << data[data.size() - 1] << "]";
    return os;
}

template <dimension_t DIM, typename SCALAR>
std::ostream& operator<<(std::ostream& os, const PhBox<DIM, SCALAR>& data) {
    os << data.min() << ":" << data.max();
    return os;
}

// Taken from boost::hash_combine
template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

}  // namespace improbable::phtree

namespace std {
template <improbable::phtree::dimension_t DIM, typename SCALAR>
struct hash<improbable::phtree::PhPointLD<DIM, SCALAR>> {
    size_t operator()(const improbable::phtree::PhPointLD<DIM, SCALAR>& x) const {
        std::size_t hash_val = 0;
        for (improbable::phtree::dimension_t i = 0; i < DIM; ++i) {
            improbable::phtree::hash_combine(hash_val, x[i]);
        }
        return hash_val;
    }
};

template <improbable::phtree::dimension_t DIM, typename SCALAR>
struct hash<improbable::phtree::PhPointHD<DIM, SCALAR>> {
    size_t operator()(const improbable::phtree::PhPointHD<DIM, SCALAR>& x) const {
        std::size_t hash_val = 0;
        for (improbable::phtree::dimension_t i = 0; i < DIM; ++i) {
            improbable::phtree::hash_combine(hash_val, x[i]);
        }
        return hash_val;
    }
};

template <improbable::phtree::dimension_t DIM, typename SCALAR>
struct hash<improbable::phtree::PhBox<DIM, SCALAR>> {
    size_t operator()(const improbable::phtree::PhBox<DIM, SCALAR>& x) const {
        std::size_t hash_val = 0;
        for (improbable::phtree::dimension_t i = 0; i < DIM; ++i) {
            improbable::phtree::hash_combine(hash_val, x.min()[i]);
            improbable::phtree::hash_combine(hash_val, x.max()[i]);
        }
        return hash_val;
    }
};
}  // namespace std
#endif  // PHTREE_COMMON_BASE_TYPES_H

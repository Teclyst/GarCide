/**
 * @file utility.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file a bunch for of utility functions used everywhere.
 * @version 0.1
 * @date 2024-07-28
 *
 * @copyright Copyright (C) 2024. Distributed under the GNU General Public
 * License, version 3.
 *
 */

/*
 * GarCide Copyright (C) 2024 Matteo Wei.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License in LICENSE for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILITY
#define UTILITY

#ifdef USE_PAR

#include <execution>

#endif

#include <iostream>
#include <ostream>
#include <regex>
#include <stdexcept>
#include <string>

namespace garcide {

/**
 * @brief Signed 8-bits integer type.
 *
 * Alias for `char`.
 *
 * We use our own types for portability.
 */
using i8 = char;

/**
 * @brief Unsigned 8-bits integer type.
 *
 * Alias for `unsigned char`.
 *
 * We use our own types for portability.
 */
using u8 = unsigned char;

/**
 * @brief Signed 16-bits integer type.
 *
 * Alias for `int`.
 *
 * We use our own types for portability.
 */
using i16 = int;

/**
 * @brief Unsigned 16-bits integer type.
 *
 * Alias for `unsigned int`.
 *
 * We use our own types for portability.
 */
using u16 = unsigned int;

/**
 * @brief Signed 32-bits integer type.
 *
 * Alias for `long`.
 *
 * We use our own types for portability.
 */
using i32 = long;

/**
 * @brief Unsigned 32-bits integer type.
 *
 * Alias for `unsigned long`.
 *
 * We use our own types for portability.
 */
using u32 = unsigned long;

/**
 * @brief Signed 64-bits integer type.
 *
 * Alias for `long long`.
 *
 * We use our own types for portability.
 */
using i64 = long long;

/**
 * @brief Unsigned 64-bits integer type.
 *
 * Alias for `unsigned long long`.
 *
 * We use our own types for portability.
 */
using u64 = unsigned long long;

/**
 * @brief Regex representing integers.
 *
 * It is a string because we often add it as a part of more complex regular
 * expressions, and there is no way to build them within `std::regex`.
 */
const std::string number_regex{"-?[1-9][0-9]*|0"};

/**
 * @brief Euclidean quotient.
 *
 * Computes the quotient in the euclidian division of `a` by `b`.
 *
 * `quot(a, b)` and `rem(a, b)` satisfy the usual relation `a == quot(a, b) * b +
 * rem(a, b)`.
 *
 * @warning Dividing by zero is undefined behaviour.
 *
 * @param a Dividend.
 * @param b Divisor.
 * @return The euclidian quotient of `a` by `b`.
 */
inline i16 quot(i16 a, i16 b) {
    i16 r = a % b;
    i16 q = a / b;
    return r >= 0 ? q : q + (b >= 0 ? -1 : 1);
};

/**
 * @brief Euclidean remainder.
 *
 * Computes the remainder in the euclidian division of `a` by `b`.
 *
 * `quot(a, b)` and `rem(a, b)` satisfy the usual relation `a == quot(a, b) * b +
 * rem(a, b)`.
 *
 * @warning Dividing by zero is undefined behaviour.
 * 
 * @param a Dividend.
 * @param b Divisor.
 * @return The euclidian remainder of `a` by `b`.
 */
inline i16 rem(i16 a, i16 b) {
    i16 r = a % b;
    return r >= 0 ? r : r + (b >= 0 ? b : -b);
};

/**
 * @brief Exception thrown in case of bad input.
 *
 * `InvalidStringError` is thrown when getting incorrect input from a function
 * that converts from a string.
 */
struct InvalidStringError {
    /**
     * @brief Source of the error.
     *
     * An error message, used to explain what went wrong.
     */
    std::string error_source;

    /**
     * @brief Construct a new `InvalidStringError` exception.
     *
     * It will hold `error_source` as
     * its error message.
     *
     * @param error_source Source of the error.
     */
    InvalidStringError(std::string error_source) : error_source(error_source) {}
};

/**
 * @brief Exception thrown when randomization is not supported.
 *
 * `NonRandomizable` is thrown when trying to randomize something for which
 * randomization is not implemented.
 */
struct NonRandomizable {};

/**
 * @brief Forward iterates applications of `f` on pairs of successive elements.
 *
 * `f` must return a boolean, and `apply_binfun` will stop at the first pair for
 * which `f` returns `true`, returning an iterator to that position.
 *
 * Linear in the sequence's length.
 *
 * @tparam ForItr A class of forward iterators.
 * @tparam BinFunc A class of binary functions of signature `bool _(&T, &T)`,
 * with `T` the type of objects `ForItr` points to.
 * @param first An iterator to the position we start iterating from.
 * @param last An iterator to the position we stop iterating at.
 * @param f A binary function of signature `bool f(&T, &T)`.
 * @return An iterator to the first position such that `f` returned `true`, or
 * to the `last` if no such position was found.
 */
template <class ForItr, class BinFunc>
inline ForItr apply_binfun(ForItr first, ForItr last, BinFunc f) {
    ForItr i, j;
    if ((i = j = first) == last)
        return first;
    while (++j != last && f(*(i++), *j))
        ;
    return j;
}

/**
 * @brief Backward iterates applications of `f` on pairs of successive elements.
 *
 * `f` must return a boolean, and `apply_binfun` will stop at the (backwardly)
 * first pair for which `f` returns `true`, returning an iterator to that
 * position.
 *
 * Linear in the sequence's length.
 *
 * @tparam BiItr A class of bidirectional iterators.
 * @tparam BinFunc A class of binary functions of signature `bool _(&T, &T)`,
 * with `T` the type of objects `BiItr` points to.
 * @param first An iterator to the position we start iterating from.
 * @param last An iterator to the position we stop iterating at.
 * @param f A binary function of signature `bool f(&T, &T)`.
 * @return An iterator to the (backwardly) first position such that `f` returned
 * `true`, or to the `last` if no such position was found.
 */
template <class BiItr, class BinFunc>
inline BiItr reverse_apply_binfun(BiItr first, BiItr last, BinFunc f) {
    BiItr i, j;
    if (first == (i = j = last))
        return first;
    --i;
    while ((j = i) != first && f(*--i, *j))
        ;
    return --j;
}

/**
 * @brief Applies `f` on pairs following the bubble sort pattern.
 *
 * _I.e._ \f[(n - 2, n - 1);
 * (n - 3, n - 2), (n - 2, n - 1); ...; (0, 1), (1, 2), ..., (n - 2, n - 1).\f]
 *
 * For each subsequence, it breaks to the next one if `f` returns `true`.
 *
 * Quadratic in the length of the sequence.
 *
 * @tparam ForItr A class of bidirectional iterators.
 * @tparam BinFun A class of binary functions of signature `bool _(&T, &T)`,
 * with `T` the type of objects `ForItr` points to.
 * @param first An iterator to the position we start bubble sorting from.
 * @param last An iterator to the position we stop bubble sorting at.
 * @param f A binary function of signature `bool f(&T, &T)`.
 */
template <class ForItr, class BinFun>
inline void bubble_sort(ForItr first, ForItr last, BinFun f) {
    ForItr i;
    if (first == (i = last))
        return;
    while (i != first)
        apply_binfun(--i, last, f);
}

/**
 * @brief A struct that represents a endline character.
 *
 * Used as an equivalent for
 * `std::endl` for `IndentedOStream`.
 *
 * You may pass an integer parameter to its constructor to specify the number of
 * lines to skip.
 */
struct EndLine {

    /**
     * @brief Number of lines to skip when this is inserted in
     * `IndentedOStream`.
     *
     * _I.e._, the number of fully white lines, as inserting an
     * `EndLine` will always result in a linebreak.
     */
    i16 lines_to_skip;

    /**
     * @brief Constructs a new `EndLine`.
     *
     * `lines_to_skip` is set to `skip`.
     *
     * @param skip What `lines_to_skip` is to be set to. Default value is 0.
     */
    EndLine(i16 skip = 0);
};

/**
 * @brief A class for output streams that keep track of indentation.
 *
 * It is not a very
 * smart implementation - the only thing it does is adding indentation after
 * seing a `Endline`. Otherwise it is just a wrapper for its `os` object.
 *
 * Because it does not actually process inputs, make sure to only ever end lines
 * with `EndLine`.
 *
 * This is not safe for multithreading.
 */
class IndentedOStream {
  private:
    /**
     * @brief Current level of indentation.
     *
     * It should be ensured that it never changes
     * after a function call.
     */
    i16 indent_level;

    /**
     * @brief The output stream `*this` is a wrapper for.
     */
    std::ostream &os;

  public:
    /**
     * @brief Constructs a new `IndentedOStream`.
     *
     * @param os The output stream `*this` wraps around.
     */
    IndentedOStream(std::ostream &os);

    /**
     * @brief Changes indentation level.
     *
     * @param new_indent_level New level of indentation.
     */
    inline void SetIndentLevel(i16 new_indent_level) {
        indent_level = new_indent_level;
    }

    /**
     * @brief Increment indentation level.
     *
     * @param indentation Variation of indentation.
     */
    inline void Indent(i16 indentation = 1) { indent_level += indentation; }

    /**
     * @brief Inserts something in the stream.
     *
     * @tparam C A class, the only requirement is that `<<<std::ostream, C>` is
     * implemented.
     * @param c What is to be inserted.
     * @return A reference to `*this`, so that `<<` may be chained.
     */
    template <class C> IndentedOStream &operator<<(const C &c) {
        os << c;
        return *this;
    }
};

/**
 * @brief Inserts an `EndLine` in the stream.
 *
 * Pushes a linebreak in `os`, then
 * flushes it and inserts the right number of spaces.
 *
 * @param el An `EndLine` value.
 * @return A reference to `*this`, so that `<<` may be chained.
 */
template <>
IndentedOStream &IndentedOStream::operator<< <EndLine>(const EndLine &el);

/**
 * @brief Indented version of `std::cout`.
 */
static IndentedOStream ind_cout(std::cout);

} // namespace garcide

#endif

/**
 * @brief Miscealanous, yet important, definitions.
 *
 * A file where we define miscealanous, yet important, functions, types and
 * classes.
 */

#ifndef UTILITY
#define UTILITY

#include <execution>
#include <iostream>
#include <ostream>
#include <regex>
#include <stdexcept>
#include <string>

namespace CGarside {

/**
 * @brief Alias for `char`.
 *
 * We use our own types for portability.
 */
typedef char sint8;
/**
 * @brief Alias for `unsigned char`.
 *
 * We use our own types for portability.
 */
typedef unsigned char uint8;
/**
 * @brief Alias for `int`.
 *
 * We use our own types for portability.
 */
typedef int sint16;
/**
 * @brief Alias for `unsigned int`.
 *
 * We use our own types for portability.
 */
typedef unsigned int uint16;
/**
 * @brief Alias for `long`.
 *
 * We use our own types for portability.
 */
typedef long sint32;
/**
 * @brief Alias for `unsigned long`.
 *
 * We use our own types for portability.
 */
typedef unsigned long uint32;
/**
 * @brief Alias for `long long`.
 *
 * We use our own types for portability.
 */
typedef long long sint64;
/**
 * @brief Alias for `unsigned long long`.
 *
 * We use our own types for portability.
 */
typedef unsigned long long uint64;

/**
 * @brief Regex representing integers.
 *
 * A regular expression that matches well-formed integers.
 * It is a string because we often add it as a part of more complex regular
 * expressions, and there is no way to build them within `std::regex`.
 */
const std::string number_regex{"-?[1-9][0-9]*|0"};

/**
 * @brief Euclidean quotient.
 *
 * Computes the quotient in the euclidian division of `a` by `b`.
 *
 * `Quot(a, b)` and `Rem(a, b)` satisfy the usual relation `a = Quot(a, b) * b +
 * Rem(a, b)`.
 *
 * @param a Dividend.
 * @param b Divisor.
 * @exception Dividing by zero is undefined behaviour.
 * @return The euclidian quotient of `a` by `b`.
 */
inline sint16 Quot(sint16 a, sint16 b) {
    sint16 r = a % b;
    sint16 q = a / b;
    return r >= 0 ? q : q + (b >= 0 ? -1 : 1);
};

/**
 * @brief Euclidean remainder.
 *
 * Computes the remainder in the euclidian division of `a` by `b`.
 *
 * `Quot(a, b)` and `Rem(a, b)` satisfy the usual relation `a = Quot(a, b) * b +
 * Rem(a, b)`.
 *
 * @param a Dividend.
 * @param b Divisor.
 * @exception Dividing by zero is undefined behaviour.
 * @return The euclidian remainder of `a` by `b`.
 */
inline sint16 Rem(sint16 a, sint16 b) {
    sint16 r = a % b;
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
     * Construct a new `InvalidStringError`, which will hold `error_source` as
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
 * Forward iterates applications of binary function `f` on pairs of successive
 * elements, between `first` and `last`.
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
 * Backward iterates applications of binary function `f` on pairs of successive
 * elements, between `first` and `last`.
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
 * Applies `f` on pairs following the bubble sort pattern (i.e. (n - 2, n - 1);
 * (n - 3, n - 2), (n - 2, n - 1); ...; (0, 1), (1, 2), ..., (n - 2, n - 1)).
 *
 * For each subsequence, it breaks to the next one if `f` returns `true`.
 *
 * Quadratic in the sequence's length.
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
 * A struct that represents a endline character. Used as an equivalent for
 * `std::endl` for `IndentedOStream`.
 *
 * You may pass an integer parameter to its constructor to specify the number of
 * lines to skip.
 */
struct EndLine {

    /**
     * @brief Number of lines to skip.
     *
     * The number of lines that should be skipped when this is inserted in
     * `IndentedOStream` (i.e., the number of fully white lines, inserting an
     * `EndLine` will always result in a linebreak).
     */
    sint16 lines_to_skip;

    /**
     * @brief Constructs a new `EndLine`.
     *
     * Constructs a new `Endline` object, whith `lines_to_skip` set to `skip`.
     *
     * @param skip What `lines_to_skip` is to be set to. Default value is 0.
     */
    EndLine(sint16 skip = 0);
};

/**
 * @brief A class for output streams that keep track of indentation.
 *
 * A class for output streams that keep track of indentation. It is not a very
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
     * The current level of indentation. You should ensure that it never changes
     * after a function call.
     */
    sint16 indent_level;

    /**
     * @brief Output stream.
     *
     * The output stream `*this` is a wrapper for.
     */
    std::ostream &os;

  public:
    /**
     * @brief Constructs a new `IndentedOStream`.
     *
     * Constructs a new `IndentedOStream`, from the output stream it wraps
     * around.
     *
     * @param os The output stream `*this` is a wrapper for.
     */
    IndentedOStream(std::ostream &os);

    /**
     * @brief Set `indent_level`.
     *
     *  Set the `indent_level` field to `new_indent_level`.
     *
     * @param indent_level New level of indentation.
     */
    inline void SetIndentLevel(sint16 new_indent_level) {
        indent_level = new_indent_level;
    }

    /**
     * @brief Increment `indent_level`.
     *
     * Increment the `indent_level` field by `indentation`.
     *
     * @param indentation Variation of indentation.
     */
    inline void Indent(sint16 indentation) { indent_level += indentation; }

    /**
     * @brief Inserts something in the stream.
     *
     * Inserts something in the stream, by inserting it in `os`.
     *
     * @tparam C A class, the only requirement is that `<<<std::ostream, C>` is
     * implemented.
     * @param c What is to be inserted.
     * @return A reference to `*this`, so that `<<` may be chained.
     */
    template <class C> IndentedOStream &operator<<(const C &c) {
        os << c;
        return *this;
    };
};

/**
 * @brief Inserts an `EndLine` in the stream.
 *
 * Inserts an `EndLine` in the stream, by going to the next line in `os`, then
 * flushing it and inserting the right amount of spaces.
 *
 * @param el An EndLine value.
 * @return A reference to `*this`, so that `<<` may be chained.
 */
template <>
IndentedOStream &IndentedOStream::operator<< <EndLine>(const EndLine &el);

/**
 * @brief Indented version of `std::cout`
 *
 * An `IndentedOStream` that wraps around `std::cout`.
 */
static IndentedOStream ind_cout(std::cout);

#ifndef USE_PAR

constexpr __pstl::execution::v1::sequenced_policy execution_policy =
    std::execution::seq;

#else

constexpr __pstl::execution::v1::parallel_policy execution_policy =
    std::execution::par;

#endif

} // namespace CGarside

#endif

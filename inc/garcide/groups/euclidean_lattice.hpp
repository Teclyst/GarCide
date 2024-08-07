/**
 * @file euclidean_lattice.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for Z ^ n.
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

#ifndef EUCLIDEAN_LATTICE
#define EUCLIDEAN_LATTICE

#include "garcide/garcide.h"

/**
 * @brief Namespace for euclidian lattices $\mathbb Z^n$.
 */
namespace garcide::euclidean_lattice {
/**
 * @brief A class for $Z^n$ canonical factors.
 *
 * A class for $Z^n$ canonical factors (elements of ${\mathbb Z / 2\mathbb
 * Z}^n$).
 */
class Underlying {

  public:
    /**
     * @brief Parameter type.
     *
     * The type for group parameters. Here it is the dimension, which is also
     * the size of the coordinates' vector, whence the use of `size_t` (type for
     * sizes).
     */
    using Parameter = size_t;

  private:
    /**
     * @brief The factor's coordinates.
     *
     * These are its coordinates in $\mathbb Z^n$ canonical basis.
     * As these are all $0$ or $1$, we use booleans, as (depending on
     * implementation) it should be more space efficient (one bit vers one
     * byte).
     */
    std::vector<bool> coordinates;

  public:
    /**
     * @brief Converts a string to a parameter.
     *
     * Converts a string to a parameter by tring to parse it as an integer.
     *
     * In case of failure, `InvalidStringError` is raised.
     *
     * @param str The string to read.
     * @return A parameter matching `str`.
     * @exception InvalidStringError Thrown in case of failure.
     */
    static Parameter parameter_of_string(const std::string &str);

    /**
     * @brief Gets the factor's dimension.
     *
     * (_I.e._ the $n$ in $\mathbb Z^n$.)
     *
     * @return The factor's dimension.
     */
    Parameter get_parameter() const;

    /**
     * @brief Access i-th coordinate.
     *
     * That way non-member functions may enjoy read read access to
     * `coordinates`.
     *
     * @param i The index that is being accessed.
     * @return The `i`-th coordinate.
     */
    inline bool at(size_t i) const { return coordinates[i]; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Its dimension will be `n`.
     *
     * `coordinates` will have length `n`, and will be filled with `false`.
     *
     * @param n The factor's `dimension`.
     */
    Underlying(Parameter n);

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, tries to extract an
     * atom to set `this` to. If it succeeds, increases `pos` so that it points
     * to just after what was extracted.
     *
     * Letting `Z = '-'? (['1' - '9'] ['0' - '9']* | '0')` be the language of
     * integers, accepted strings are those represented by regular expression
     * `('e' '_'?)? Z`, under the additional hypothesis that the integer they
     * represent is in [`0`, `dimension`[, ignoring whitespaces.
     *
     * "e_i" (with optional "e_") stands for base vector e_i.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches `Z`, or if there is one, if the corresponding
     * integer does not belong to [`0`, `get_parameter()`[.
     */
    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Height of the lattice.
     *
     * (_I.e._ `delta`'s length as a word in the generators, here
     * `dimension`.)
     *
     * @return The height of the lattice.
     */
    inline sint16 lattice_height() const { return int(get_parameter()); }

    /**
     * @brief Prints internal data in `os`.
     *
     * Prints the factor's `coordinates` in `os`, typically for
     * debugging purposes.
     *
     * @param os The output stream it is printed in.
     */
    void debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor in `os`.
     *
     * It is printed as a product (multiplicative convention) of
     * base vectors.
     *
     * @param os The output stream it is printed in.
     */
    void print(IndentedOStream &os) const;

    /**
     * @brief Sets the factor to the identity.
     *
     * (_I.e._ sets all coordinates to `false`.)
     */
    void identity();

    /**
     * @brief Sets the factor to delta.
     *
     * (_I.e._ sets all coordinates to `true`.)
     */
    void delta();

    /**
     * @brief Computes meets.
     *
     * Computes the (left, although that does not matter here as Z ^ n is
     * abelian) meet of `*this` and `b` (i.e., coordinates-wise `&&`).
     *
     * @param b Second argument.
     * @return The meet of `*this` and `b`.
     */
    Underlying left_meet(const Underlying &b) const;

    /**
     * @brief Computes meets.
     *
     * Computes the (right, although that does not matter here as Z ^ n is
     * abelian) meet of `*this` and `b` (i.e., coordinates-wise `&&`).
     *
     * @param b Second argument.
     * @return The meet of `*this` and `b`.
     */
    inline Underlying right_meet(const Underlying &b) const {
        return left_meet(b);
    };

    /**
     * @brief Equality check.
     *
     * Compares `*this` and `b`, returning `true` if they are equal (i.e.
     * coordinate-wise equal).
     *
     * @param b Second argument.
     * @return If `*this` and `b` are equal.
     */
    inline bool compare(const Underlying &b) const {
        return coordinates == b.coordinates;
    };

    /**
     * @brief product computations.
     *
     * Computes the product of `*this` and `b` (i.e. coordinate-wise xor), under
     * the assumption that it lies below delta (this is not actually checked).
     *
     * @param b Second argument.
     * @return The product of `*this` and `b`.
     */
    Underlying product(const Underlying &b) const;

    /**
     * @brief Complement computations.
     *
     * Computes the complement of `*this` to `b` (i.e. the factor `c` such
     * that `ac=b`), under the assumption that `*this` is smaller than `b`. As Z
     * ^ n is abelian, left and right variants are in fact the same.
     *
     * In this case, this is actually the same (well, as an operation on bit
     * vectors) as product.
     *
     * @param b Second argument.
     * @return The complement of `*this` to `b`.
     */
    inline Underlying left_complement(const Underlying &b) const {
        return product(b);
    };

    /**
     * @brief Complement computations.
     *
     * Computes the complement of `*this` to `b` (i.e. the factor `c` such
     * that `ac = b`), under the assumption that `*this` is smaller than `b`. As
     * Z ^ n is abelian, left and right variants are in fact the same.
     *
     * In this case, this is actually the same (well, as an operation on bit
     * vectors) as product.
     *
     * @param b Second argument.
     * @return The complement of `*this` to `b`.
     */
    inline Underlying right_complement(const Underlying &b) const {
        return product(b);
    };

    /**
     * @brief Sets `*this` to a random factor.
     *
     * Sets `*this` to a random factor, following an uniform distribution other
     * factors.
     */
    void randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (i.e. base vectors).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by delta ^ k.
     *
     * Conjugates `*this` by delta ^ k (actually, does nothing as Z ^ n is
     * abelian).
     *
     * @param k The exponent.
     */
    inline void delta_conjugate_mut(__attribute__((unused)) sint16 k) {};

    /**
     * @brief Hashes the factor.
     *
     * Hashes the factor. Done by interpreting it as a polynomial and evaluating
     * it in 2 (yields a bijection between the set of factors, and [0, 2 ^ n[).
     *
     * @return The hash.
     */
    size_t hash() const;
};

/**
 * @brief Class for canonical factors.
 *
 * Class for canonical factors, will all corresponding methods.
 */
typedef FactorTemplate<Underlying> Factor;

/**
 * @brief Class for elements of Z ^ n.
 *
 * Class for elements of Z ^ n, as a Garside group.
 */
typedef BraidTemplate<Factor> Braid;

} // namespace garcide::euclidean_lattice

#endif
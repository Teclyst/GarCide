/**
 * @file euclidean_lattice.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for euclidian lattices \f$\mathbb Z^n\f$.
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

#include "garcide/garcide.hpp"

/**
 * @brief Namespace for euclidian lattices \f$\mathbb Z^n\f$.
 */
namespace garcide::euclidean_lattice {

/**
 * @brief A class for \f$\mathbb Z^n\f$ canonical factors.
 *
 * They are represented by elements of \f${\mathbb Z
 * / 2\mathbb Z}^n\f$, as bitvectors.
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
     * These are its coordinates in the canonical basis of \f$\mathbb Z^n\f$.
     * As these are all \f$0\f$ or \f$1\f$, we use booleans, as (depending on
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
     * (_I.e._ the \f$n\f$ in \f$\mathbb Z^n\f$.)
     *
     * @return The factor's dimension.
     */
    Parameter get_parameter() const;

    /**
     * @brief Access i-th coordinate (read-only).
     *
     * @param i The index that is being accessed.
     * @return The `i`-th coordinate.
     */
    inline bool at(size_t i) const { return coordinates[i]; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Its dimension will be `n`, and it will be initialized as a `n`-length
     * vector filled with `false`.
     *
     * @param n The dimension.
     */
    Underlying(Parameter n);

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, tries to extract an
     * atom to set `this` to. If it succeeds, increases `pos` so that it points
     * to just after what was extracted.
     *
     * Letting \f$Z = \texttt{-}? ([\texttt{1} - \texttt{9}] [\texttt{0} -
     * \texttt{9}]^* \mid \texttt{0})\f$ be a regular expression representing
     * integers, accepted strings are those represented by regular expression
     * \f$(\texttt{e} \texttt{_}?)? Z\mid \texttt{D}\f$, under the additional hypothesis that the integer they
     * represent is in \f$[0, n[\f$, where \f$n\f$ is the dimension, and ignoring
     * whitespaces.
     *
     * \f$\texttt{e_}i\f$, \f$\texttt{e}i\f$ and \f$i\f$ all stand for base vector \f$e_i\f$.
     *
     * \f$\texttt{D}\f$  represents \f$\Delta\f$.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches the expression, or if there is one, if the corresponding
     * integer does not belong to \f$[0, n[\f$.
     */
    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Height of the lattice.
     *
     * (_I.e._ the length of \f$\Delta\f$ as a word in the generators, here
     * the dimension.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const { return int(get_parameter()); }

    /**
     * @brief Prints internal data in `os`.
     *
     * Prints `private` member `coordinates`, typically for debugging.
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
     *
     * Linear in the dimension.
     */
    void identity();

    /**
     * @brief Sets the factor to \f$\Delta\f$.
     *
     * (_I.e._ sets all coordinates to `true`.)
     *
     * Linear in the dimension.
     */
    void delta();

    /**
     * @brief Computes meets.
     *
     * Computes the (left, although that does not matter here as \f$\mathbb Z ^
     * n\f$ is abelian) meet of `*this` and `b` (_i.e._, coordinates-wise `&&`).
     *
     * Linear in the dimension.
     *
     * @param b Second operand.
     * @return The meet of `*this` and `b`.
     */
    Underlying left_meet(const Underlying &b) const;

    /**
     * @brief Computes meets.
     *
     * Computes the (right, although that does not matter here as \f$\mathbb Z ^
     * n\f$ is abelian) meet of `*this` and `b` (_i.e._, coordinates-wise `&&`).
     *
     * Linear in the dimension.
     *
     * @param b Second operand.
     * @return The meet of `*this` and `b`.
     */
    inline Underlying right_meet(const Underlying &b) const {
        return left_meet(b);
    };

    /**
     * @brief Equality check.
     *
     * Compares `*this` and `b`, returning `true` if they are equal (_i.e._
     * coordinate-wise equal).
     *
     * Linear in the dimension.
     *
     * @param b Second operand.
     * @return If `*this` and `b` are equal.
     */
    inline bool compare(const Underlying &b) const {
        return coordinates == b.coordinates;
    };

    /**
     * @brief Product computations.
     *
     * Computes the product of `*this` and `b` (_i.e._ coordinate-wise xor),
     * under the assumption that it is a factor.
     *
     * Linear in the dimension.
     *
     * @param b Second operand.
     * @return The product of `*this` and `b`.
     */
    Underlying product(const Underlying &b) const;

    /**
     * @brief Complement computations.
     *
     * Computes the complement of `*this` to `b` under the assumption that
     * `*this` is smaller than `b`. As \f$\mathbb Z ^ n\f$ is abelian, left and
     * right variants are in fact the same.
     *
     * In this case, this is actually the same (as an operation on bit
     * vectors) as product.
     *
     * Linear in the dimension.
     *
     * @param b Second operand.
     * @return The complement of `*this` to `b`.
     */
    inline Underlying left_complement(const Underlying &b) const {
        return product(b);
    };

    /**
     * @brief Complement computations.
     *
     * Computes the complement of `*this` to `b` under the assumption that
     * `*this` is smaller than `b`. As \f$\mathbb Z ^ n\f$ is abelian, left and
     * right variants are in fact the same.
     *
     * In this case, this is actually the same (as an operation on bit
     * vectors) as product.
     *
     * Linear in the dimension.
     *
     * @param b Second operand.
     * @return The complement of `*this` to `b`.
     */
    inline Underlying right_complement(const Underlying &b) const {
        return product(b);
    };

    /**
     * @brief Sets `*this` to a random factor.
     *
     * It is chosen uniformly other \f$(\mathbb Z/2\mathbb Z)^n\f$.
     *
     * Linear in the dimension.
     */
    void randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (_i.e._ base vectors).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by \f$\Delta^k\f$.
     *
     * Actually, this does nothing as \f$\mathbb Z^n\f$ is abelian.
     *
     * Constant time.
     *
     * @param k The exponent.
     */
    inline void delta_conjugate_mut(__attribute__((unused)) i16 k) {};

    /**
     * @brief Hashes the factor.
     *
     * Done by interpreting it as a polynomial and evaluating
     * it in \f$2\f$ (yields a bijection between the set of factors and \f$[0,
     * 2 ^ n[\f$).
     *
     * Linear in the dimension.
     *
     * @return The hash.
     */
    size_t hash() const;
};

/**
 * @brief Class for \f$\mathbb Z^n\f$ canonical factors, as a Garside group.
 */
using Factor = FactorTemplate<Underlying>;

/**
 * @brief Class for \f$\mathbb Z ^ n\f$ elements, as a Garside group.
 */
using Braid = BraidTemplate<Factor>;

} // namespace garcide::euclidean_lattice

#endif
/**
 * @file artin.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for standard braid groups (classic Garside structure).
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

#ifndef ARTIN
#define ARTIN

#include "garcide/ultra_summit.h"

namespace garcide {

/**
 * @brief Namespace for standard braid groups, classic Garside structure.
 *
 * Contains the usual `Underlying` class, as well as functions to compute
 * Thurston types.
 */
namespace artin {

/**
 * @brief A class for classic Garside structure braid group canonical factors.
 *
 * They are represented by permutations, stored as inverse permutation tables.
 */
class Underlying {

  public:
    /**
     * @brief Parameter type.
     *
     * Here its elements represent positive integers.
     */
    using Parameter = i16;

  private:
    /**
     * @brief Inverse table of the permutation that corresponds to the factor.
     *
     * Given a permutation \f$\sigma\in\mathfrak S_n\f$, the inverse permutation
     * table of \f$\sigma\f$ is the table $\f\mathrm T_\sigma\f$ such that
     * $\forall i\in[\![1,n]\!],\mathrm T_\sigma\f[\sigma(i)]=i\f$.
     *
     * Indexes start at \f$1\f$.
     */
    std::vector<i16> permutation_table;

  public:
    /**
     * @brief Maximum number of strands.
     *
     * The greatest number of strands that may be used for braids.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const i16 MAX_NUMBER_OF_STRANDS = 256;

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
     * @brief Gets the number of strands.
     *
     * That is to say, the \f$n\f$ in \f$B_n\f$.
     *
     * @return The number of strands.
     */
    inline Parameter get_parameter() const {
        return int(permutation_table.size() - 1);
    }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `n` as its
     * number of strands.
     *
     * Its permutation table will have length `n`, and will be filled with
     * zeros (thus this is not a valid factor). It should be initialized it with
     * `identity()`, `delta()`, or another similar member.
     *
     * @param n The number of strands of the factor (also the length of
     * its permutation table).
     */
    Underlying(Parameter n);

    /**
     * @brief Access the `i`-th element of the permutation table (read-only).
     *
     * `i` should be between `1` and the number of strands.
     *
     * @param i The index that is being accessed.
     * @return The `i`-th element of the permutation table.
     */
    inline i16 at(size_t i) const { return permutation_table[i]; }

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting \f$Z = \texttt{-}? ([\texttt{1} - \texttt{9}] [\texttt{0} -
     * \texttt{9}]^* \mid \texttt{0})\f$ be a regular expression representing
     * integers, accepted strings are those represented by regular expression
     * \f$(\texttt{s} \texttt{_}?)? Z\mid \texttt{D}\f$, under the additional
     * hypothesis that the integer they represent is in \f$[1, n[\f$, where
     * \f$n\f$ is the number of strands, and ignoring whitespaces.
     *
     * \f$\texttt{s_}i\f$, \f$\texttt{s}i\f$ and \f$i\f$ all stand for
     * Artin generator \f$\sigma_i\f$.
     *
     * \f$\texttt{D}\f$ represents \f$\Delta\f$.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches the expression, or if there is one, if the
     * corresponding integer does not belong to \f$[1,n[\f$.
     */
    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Height of the lattice.
     *
     * (_I.e._, letting `n` be the number of of strands, the length of
     * \f$\Delta_n\f$ as a word in the generators, which is
     * \f$\frac{n(n-1)}2\f$.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const {
        Parameter n = get_parameter();
        return n * (n - 1) / 2;
    }

    /**
     * @brief Prints internal representation in `os`.
     *
     * Prints private member `permutation_table`, typically for debugging.
     *
     * @param os The output stream it is printed in.
     */
    void debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor to `os`.
     *
     * It is printed as a product in the Artin generators.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    /**
     * @brief Sets the factor to the identity.
     *
     * (_I.e._ sets the permutation table to the one of the identity
     * permutation.)
     *
     * Linear in the number of strands.
     */
    void identity();

    /**
     * @brief Sets the factor to the Garside element.
     *
     * (_I.e._ sets the permutation table to the one of \f$\Delta_n\f$,
     * where \f$n\f$ is the number of strands. This is the table \f$\mathrm
     * T_{\Delta_n}\f$ such that \f$\forall
     * i\in[\![1,n]\!],T_{\Delta_n}[i]=n-i+1\f$.)
     *
     * Linear in the number of strands.
     */
    void delta();

    /**
     * @brief Computes the left meet of `*this` and `b`.
     *
     * Uses a recursive algorithm by Thurston, that is presented in
     * Cha, Ko, Lee, Han, Cheon, _An Efficient Implementation of Braid Groups_,
     * 2001.
     *
     * Quasilinear in the number of strands.
     *
     * @param b Second operand.
     * @return The left meet of `*this` and `b`.
     */
    Underlying left_meet(const Underlying &b) const;

    /**
     * @brief Computes the right meet of `*this` and `b`.
     *
     * Uses a recursive algorithm by Thurston, that is presented in
     * Cha, Ko, Lee, Han, Cheon, _An Efficient Implementation of Braid Groups_,
     * 2001.
     *
     * Quasilinear in the number of strands.
     *
     * @param b Second operand.
     * @return The right meet of `*this` and `b`.
     */
    Underlying right_meet(const Underlying &b) const;

    /**
     * @brief Equality check.
     *
     * Compares `*this` and `b`, returning `true` if they are equal (_i.e._
     * they represent the same permutation).
     *
     * Linear in the number of strands.
     *
     * @param b Second operand.
     * @return If `*this` and `b` are equal.
     */
    inline bool compare(const Underlying &b) const {
        return permutation_table == b.permutation_table;
    }

    /**
     * @brief Product computations.
     *
     * Computes the product of `*this` and `b` (_i.e._ composition of
     * permutations), under the assumption that it is a factor.
     *
     * Linear in the number of strands.
     *
     * @param b Second (right) operand.
     * @return The product of `*this` and `b`.
     */
    Underlying product(const Underlying &b) const;

    /**
     * @brief Left complement computations.
     *
     * Computes the left complement of `*this` to `b`, under the assumption that
     * `*this` left-divides `b`.
     *
     * Linear in the number of strands.
     *
     * @param b Second operand.
     * @return The left complement of `*this` to `b`.
     */
    Underlying left_complement(const Underlying &b) const;

    /**
     * @brief Right complement computations.
     *
     * Computes the right complement of `*this` to `b`, under the assumption
     * that `*this` right-divides `b`.
     *
     * Linear in the number of strands.
     *
     * @param b Second operand.
     * @return The right complement of `*this` to `b`.
     */
    Underlying right_complement(const Underlying &b) const;

    /**
     * @brief Sets `*this` to a random factor.
     *
     * It is chosen uniformly other \f$\mathfrak S_n\f$, where \f$n\f$ is the
     * dimension, using Knuth's shuffle algorithm.
     *
     * Linear in the number of strands.
     */
    void randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (_i.e._ the Artin generators).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by \f$\Delta_n^k\f$.
     *
     * Here \f$n\f$ denotes the number of strands.
     *
     * Linear in the number of strands (but does nothing if `k` is even).
     *
     * @param k The exponent.
     */
    void delta_conjugate_mut(i16 k);

    /**
     * @brief Hashes the factor.
     *
     * Linear in the number of strands.
     *
     * @return The hash.
     */
    size_t hash() const;

    /**
     * @brief Computes the tableau associated with a factor.
     *
     * Computes the tableau associated with `this` and store it in `tab`.
     *
     * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses'
     * code.
     *
     * @param tab A matrix where the tableau is to be stored.
     */
    void tableau(i16 **&tab) const;

  private:
    /**
     * @brief Computes the factor associated with the inverse of the permutation
     * of this factor.
     *
     * Linear in the number of strands.
     *
     * @return The factor associated with the inverse of the permutation of
     * `*this`.
     */
    Underlying inverse() const;

    // Subroutine called by left_meet() and right_meet().
    static void MeetSub(const i16 *a, const i16 *b, i16 *r, i16 s,
                        i16 t);
};

/**
 * @brief Class for classic Garside structure braid groups canonical factors.
 */
using Factor = FactorTemplate<Underlying>;

/**
 * @brief Class for classic Garside structure braid groups elements.
 */
using Braid = BraidTemplate<Factor>;

/**
 * @brief `enum` for Thurston types.
 *
 * An `enum` whose elements represent the three Thurston types.
 */
enum class ThurstonType { Periodic, Reducible, PseudoAsonov };

/**
 * @brief Determines if a braid preserves a family of circles.
 *
 * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses' code.
 *
 * @param b The braid to be tested.
 * @return Whether `b` preserves a family of circles.
 */
bool preserves_circles(const Braid &b);

/**
 * @brief Computes the Thurston type of a braid whose USS was already computed.
 *
 * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses' code.
 *
 * @param b The braid whose Thurston type is to be computed.
 * @param uss The ultra summit set of `b`.
 * @return The Thurston type of `b`.
 */
ThurstonType thurston_type(const Braid &b,
                           const ultra_summit::UltraSummitSet<Braid> &uss);

/**
 * @brief Computes the Thurston type of a braid.
 *
 * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses' code.
 *
 * @param b The braid whose Thurston type is to be computed.
 * @return the Thurston type of `b`.
 */
ThurstonType thurston_type(const Braid &b);

} // namespace artin

/**
 * @brief Inserts a Thurston type in the output stream.
 *
 * @param type The Thurston type to be inserted.
 * @return A reference to `*this`, so that `<<` may be chained.
 */
template <>
IndentedOStream &IndentedOStream::operator<< <artin::ThurstonType>(
    const artin::ThurstonType &type);

} // namespace garcide

#endif
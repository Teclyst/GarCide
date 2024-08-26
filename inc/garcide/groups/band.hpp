/**
 * @file band.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for standard braid groups (dual Garside structure).
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

#ifndef BAND
#define BAND

#include "garcide/garcide.h"

#ifdef USE_CLN

#include <cln/integer.h>

#endif

/**
 * @brief Namespace for standard braid groups, dual Garside structure.
 */
namespace garcide::band {

/**
 * @brief A class for dual Garside structure braid group canonical factors.
 *
 * They are represented by permutations, stored as inverse permutation tables.
 *
 * To a factor one can also associate a partition of \f$[\![1,n]\!]\f$, where
 * \f$n\f$ is the number of strands. These are represented by integer arrays,
 * with the following convention: given an array \f$\mathrm A\f$ representing a
 * partition, then for all \f$i\in[\![1,n]\!]\f$, \f$\mathrm A[i]\f$ is the
 * minimum of the cell \f$i\f$ belongs to.
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
    static const Parameter MAX_NUMBER_OF_STRANDS = 256;

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
     * @brief Height of the lattice.
     *
     * (_I.e._, letting `n` be the number of of strands, the length of
     * \f$\delta_n\f$ as a word in the generators, which is
     * \f$n-1\f$.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const { return get_parameter() - 1; }

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
    Underlying(i16 n);

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
     * Reads the string `str`, starting at position `pos`, tries to extract an
     * atom to set `this` to. If it succeeds, increases `pos` so that it points
     * to just after what was extracted.
     *
     * Letting \f$Z = \texttt{-}? ([\texttt{1} - \texttt{9}] [\texttt{0} -
     * \texttt{9}]^* \mid \texttt{0})\f$ be a regular expression representing
     * integers, accepted strings are those represented by regular expression
     * \f$(\texttt{a} \texttt{_}?)?\texttt{(}Z \texttt{,}? Z\texttt{)}\mid
     * \texttt{D}\f$, under the additional assumptions that the integers they
     * represent are in \f$[1, n]\f$, where \f$n\f$ is the number of strands,
     * and not equal. Whitespaces are ignored.
     *
     * \f$\texttt{a_(}i\texttt{,}j\texttt)\f$,
     * \f$\texttt{a_(}i\texttt{ }j\texttt)\f$,
     * \f$\texttt{a(}i\texttt{,}j\texttt)\f$,
     * \f$\texttt{a(}i\texttt{ }j\texttt)\f$,
     * \f$\texttt{(}i\texttt{,}j\texttt)\f$ and
     * \f$\texttt{(}i\texttt{ }j\texttt)\f$
     * all stand for Birman-Ko-Lee generator \f$a_{i,j}\f$.
     *
     * \f$\texttt{D}\f$ represents \f$\delta_n\f$.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches the expression, or if there is one, if
     * either integer does not belong to \f$[1,n]\f$, or if they are equal.
     */
    void of_string(const std::string &str, size_t &pos);

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
     * It is printed as a product in the Birman-Ko-Lee generators, in a way that
     * spells the disjunct cycle decomposition.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    /**
     * @brief Computes the partition associated with the factor.
     *
     * It is then written in `i16` array `x`.
     *
     * A partition is represented by an integer array,
     * with the following convention: given an array \f$\mathrm A\f$
     * representing a partition, then for all \f$i\in[\![1,n]\!]\f$, \f$\mathrm
     * A[i]\f$ is the minimum of the cell \f$i\f$ belongs to.
     *
     * Linear in the number of strands.
     * 
     * @param x A `i16` array.
     */
    void assign_partition(i16 *x) const;

    /**
     * @brief Sets the factor to one that corresponds to the partition stored in
     * `x`.
     *
     * A partition is represented by an integer array,
     * with the following convention: given an array \f$\mathrm A\f$
     * representing a partition, then for all \f$i\in[\![1,n]\!]\f$, \f$\mathrm
     * A[i]\f$ is the minimum of the cell \f$i\f$ belongs to.
     *
     * Linear in the number of strands.
     * 
     * @param x A `i16` array.
     */
    void of_partition(const i16 *x);

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
     * (_I.e._ sets the permutation table to the one of \f$\delta_n\f$,
     * where \f$n\f$ is the number of strands. This is the table \f$\mathrm
     * T_{\delta_n}\f$ such that \f$\forall
     * i\in[\![1,n]\!],T_{\delta_n}[i]\equiv i+1\ [n]\f$.)
     *
     * Linear in the number of strands.
     */
    void delta();

    /**
     * @brief Computes the meet of `*this` and `b`.
     *
     * This is done by getting the partitions associated with them, then
     * computing their meet.
     *
     * For the dual structure, the left and right meets are equal.
     *
     * Linear in the number of strands, although it uses a `thread_local` square
     * matrix of dimension `MAX_BRAID_INDEX + 1`.
     *
     * @param b Second operand.
     * @return The meet of `*this` and `b`.
     */
    Underlying left_meet(const Underlying &b) const;

    /**
     * @brief Computes the meet of `*this` and `b`.
     *
     * This is done by getting the partitions associated with them, then
     * computing their meet.
     *
     * For the dual structure, the left and right meets are equal.
     *
     * Linear in the number of strands, although it uses a `thread_local` square
     * matrix of dimension `MAX_BRAID_INDEX + 1`.
     *
     * @param b Second operand.
     * @return The meet of `*this` and `b`.
     */
    inline Underlying right_meet(const Underlying &b) const {
        return left_meet(b);
    }

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
     * `*this` divides `b`.
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
     * that `*this` divides `b`.
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
     * This is done by polling a random ballot sequence, which are in bijection
     * with noncrossing partitions (see Cha, Ko, Lee, Han, Cheon, _An Efficient
     * Implementation of Braid Groups_, 2001).
     *
     * However, polling a random ballot sequence is done essentially though
     * enumerating them. Therefore we need bigints, and the original _CBraid_
     * code uses the _CLN_ library, that is not easy to set up with _CMake_. It
     * is planned to switch from _CLN_ to a more specific bigints library, but
     * for the time being the function will just throw `NonRandomizable`.
     *
     * It should end up being quasilinear in the number of strands.
     *
     * @exception NonRandomizable Thrown, unless by some miracle you are
     * compiling the project with _CLN_.
     */
    void randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (_i.e._ the Birman-Ko-Lee generators).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by \f$\Delta_n^k\f$.
     *
     * Here \f$n\f$ denotes the number of strands.
     *
     * Linear in the number of strands (in particular, does not depend on `k`).
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
     * @brief Sets the factor to the one associated with a given ballot
     * sequence.
     *
     * @param s The ballot sequence.
     */
    void of_ballot_sequence(const i8 *s);

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
};

/**
 * @brief Class for dual Garside structure braid groups canonical factors.
 */
using Factor = FactorTemplate<Underlying>;

/**
 * @brief Class for dual Garside structure braid groups elements.
 */
using Braid = BraidTemplate<Factor>;

#ifdef USE_CLN

/**
 * @brief Sets `s` to the `k`-th ballot sequence of length `2 * n`.
 *
 * The order (and the algorithm) are described in Cha, Ko, Lee, Han, Cheon, _An
 * Efficient Implementation of Braid Groups_, 2001
 *
 * @param n Length of the ballot sequence, divided by \f$2\f$.
 * @param k Rank of the ballot sequence.
 * @param s Integer array to be set to the ballot sequence.
 */
void ballot_sequence(i16 n, cln::cl_I k, i8 *s);

/**
 * @brief Returns the `n`-th Catalan number.
 *
 * These are pre-computed, so this is a constant-time operation.
 *
 * @param n Rank of the Catalan number.
 * @return The `n`-th Catalan number.
 */
const cln::cl_I &get_catalan_number(i16 n);

#endif

} // namespace garcide::band

#endif
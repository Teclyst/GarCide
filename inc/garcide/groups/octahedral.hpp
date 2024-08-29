/**
 * @file octahedral.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for \f$\mathbf B\f$-series Artin groups (dual Garside
 * structure).
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

#ifndef OCTAHEDRAL
#define OCTAHEDRAL

#include "garcide/garcide.h"

/**
 * @brief Namespace for \f$\mathbf B\f$-series Artin groups, dual Garside
 * structure.
 */
namespace garcide::octahedral {

/**
 * @brief A class for dual Garside structure \f$\mathbf B\f$-series Artin groups
 * canonical factors.
 *
 * They are represented by permutations (of \f$[\![1,2n]\!]\f$), stored as
 * inverse permutation tables.
 *
 * To a factor one can also associate a partition of \f$[\![1,2n]\!]\f$, where
 * \f$n\f$ is the number of strands. These are represented by integer arrays,
 * with the following convention: given an array \f$\mathrm A\f$ representing a
 * partition, then for all \f$i\in[\![1,2n]\!]\f$, \f$\mathrm A[i]\f$ is the
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
     * Given a permutation \f$\sigma\in\mathfrak S_{2n}\f$, the inverse
     * permutation table of \f$\sigma\f$ is the table $\f\mathrm T_\sigma\f$
     * such that
     * $\forall i\in[\![1,2n]\!],\mathrm T_\sigma\f[\sigma(i)]=i\f$.
     *
     * Indexes start at \f$1\f$.
     */
    std::vector<i16> permutation_table;

  public:
    /**
     * @brief Maximum group parameter.
     *
     * The greatest group parameter (_i.e._, the \f$n\f$ in \f$\mathrm
     * A_n(\mathbf B)\f$) that may be used for factors.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const i16 MAX_PARAMETER = 256;

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
     * @brief Gets the group parameter.
     *
     * That is to say, the \f$n\f$ in \f$\mathrm A_n(\mathbf B)\f$.
     *
     * @return The group parameter.
     */
    inline i16 get_parameter() const { return permutation_table.size() / 2; }

    /**
     * @brief Height of the lattice.
     *
     * (_I.e._, letting `n` be the parameter, the length of
     * \f$\Delta\f$ as a word in the generators, which is
     * \f$n\f$.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const { return get_parameter(); }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `n` as its
     * parameter.
     *
     * Its permutation table will have length `2 * n`, and will be filled with
     * zeros (thus this is not a valid factor). It should be initialized it with
     * `identity()`, `delta()`, or another similar member.
     *
     * @param n The parameter of the factor (also the length of
     * its permutation table, divided by \f$2\f$).
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
     * \f$(\texttt{s} \texttt{_}?)?\texttt{(}Z \texttt{,}? Z\texttt{)}\mid
     * (\texttt{l} \texttt{_}?)? Z \mid \texttt{D}\f$, under the additional
     * assumption that the integers they represent are not equal mod \f$n\f$,
     * where \f$n\f$ is the parameter. Whitespaces are ignored.
     *
     * \f$\texttt{s_(}i\texttt{,}j\texttt)\f$,
     * \f$\texttt{s_(}i\texttt{ }j\texttt)\f$,
     * \f$\texttt{s(}i\texttt{,}j\texttt)\f$,
     * \f$\texttt{s(}i\texttt{ }j\texttt)\f$,
     * \f$\texttt{(}i\texttt{,}j\texttt)\f$ and
     * \f$\texttt{(}i\texttt{ }j\texttt)\f$
     * all stand for the Bessis **short** generator that swaps \f$i\f$ and \f$j\f$,
     * mod \f$2n\f$.
     *
     * \f$\texttt{l_}i\f$,
     * \f$\texttt{l}i\f$ and
     * \f$i\f$ represent the Bessis **long** generator that swaps \f$i\f$ and
     * \f$i+n\f$, mod \f$2n\f$.
     *
     * See Bessis, _The Dual Braid Monoid_, 2001, [arXiv:math/0101158
     * [math.GR]](https://arxiv.org/abs/math/0101158).
     *
     * \f$\texttt{D}\f$ represents \f$\Delta\f$.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches the expression, or if there is one,
     * it matches the first part, and both integers are equal mod \f$n\f$.
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
     * It is printed as a product in the Bessis short and long generators, in a
     * way that spells the disjunct cycle decomposition.
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
     * representing a partition, then for all \f$i\in[\![1,2n]\!]\f$, \f$\mathrm
     * A[i]\f$ is the minimum of the cell \f$i\f$ belongs to.
     *
     * Linear in the parameter.
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
     * representing a partition, then for all \f$i\in[\![1,2n]\!]\f$, \f$\mathrm
     * A[i]\f$ is the minimum of the cell \f$i\f$ belongs to.
     *
     * Linear in the parameter.
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
     * Linear in the parameter.
     */
    void identity();

    /**
     * @brief Sets the factor to the Garside element.
     *
     * (_I.e._ sets the permutation table to the one of \f$\Delta\f$. This is
     * the table \f$\mathrm T_{\Delta}\f$ such that \f$\forall
     * i\in[\![1,2n]\!],T_{\Delta}[i]\equiv i+1\ [2n]\f$.)
     *
     * Linear in the parameter.
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
     * Linear in the parameter, although it uses a `thread_local` square
     * matrix of dimension `2 * MAX_PARAMETER + 1`.
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
     * Linear in the parameter, although it uses a `thread_local` square
     * matrix of dimension `2 * MAX_PARAMETER + 1`.
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
     * Linear in the parameter.
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
     * Linear in the parameter.
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
     * Linear in the parameter.
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
     * Linear in the parameter.
     *
     * @param b Second operand.
     * @return The right complement of `*this` to `b`.
     */
    Underlying right_complement(const Underlying &b) const;

    /**
     * @brief Sets `*this` to a random factor.
     *
     * Currently always throws `NonRandomizable`.
     *
     * @exception NonRandomizable Thrown.
     */
    inline void randomize() { throw NonRandomizable(); }

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (_i.e._ Bessis short and long generators).
     *
     * See Bessis, _The Dual Braid Monoid_, 2001, [arXiv:math/0101158
     * [math.GR]](https://arxiv.org/abs/math/0101158).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by \f$\Delta^k\f$.
     *
     * Linear in the parameter (in particular, does not depend on `k`).
     *
     * @param k The exponent.
     */
    void delta_conjugate_mut(i16 k);

    /**
     * @brief Hashes the factor.
     *
     * Linear in the parameter.
     *
     * @return The hash.
     */
    std::size_t hash() const;

  private:
    /**
     * @brief Computes the factor associated with the inverse of the permutation
     * of this factor.
     *
     * Linear in the parameter.
     *
     * @return The factor associated with the inverse of the permutation of
     * `*this`.
     */
    Underlying inverse() const;
};

/**
 * @brief Class for dual Garside structure \f$\mathbf B\f$-series Artin groups
 * canonical factors.
 */
typedef FactorTemplate<Underlying> Factor;

/**
 * @brief Class for dual Garside structure \f$\mathbf B\f$-series Artin groups
 * elements.
 */
typedef BraidTemplate<Factor> Braid;

} // namespace garcide::octahedral

#endif
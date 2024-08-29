/**
 * @file dual_complex.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for \f$\mathrm B(e,e,n+1)\f$ groups (dual Garside
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

#include "garcide/garcide.h"

namespace garcide {

/**
 * @brief Namespace for complex braid groups \f$\mathrm B(e,e,n+1)\f$, dual
 * Garside structure.
 */
namespace dual_complex {

/**
 * @brief Class for complex braid groups parameters.
 *
 * An instance with members `e` and `n` set to \f$e\f$ and \f$n\f$ respectively
 * is the parameter of \f$\mathrm B(e,e,n+1)\f$.
 */
struct EENParameter {
    /**
     * @brief \f$e\f$ field.
     *
     * It is the `e` such that coefficients are in \f$\mathbb U_e\f$.
     */
    i16 e;

    /**
     * @brief \f$n\f$ field.
     *
     * Factors are (conceptually) square matrix of dimension `n+1`.
     */
    i16 n;

    /**
     * @brief Construct a new `EENParameter`.
     *
     * @param e Its `e` member.
     * @param n Its `n` member.
     */
    inline EENParameter(i16 e, i16 n) : e(e), n(n) {}

    /**
     * @brief Equality check.
     *
     * @param p Second operand.
     * @return If `*this` and `p` are equal.
     */
    inline bool compare(const EENParameter &p) const {
        return ((e == p.e) && (n == p.n));
    }

    /**
     * @brief Equality check.
     *
     * Syntactic sugar for `compare()`.
     *
     * @param p Second operand.
     * @return If `*this` and `p` are equal.
     */
    inline bool operator==(const EENParameter &p) const { return compare(p); }

    /**
     * @brief Unequality check.
     *
     * @param p Second operand.
     * @return If `*this` and `p` are not equal.
     */
    inline bool operator!=(const EENParameter &p) const { return !compare(p); }

    /**
     * @brief Prints the parameter to output stream `os`.
     *
     * @param os The output stream it is printed in.
     */
    void print(IndentedOStream &os) const;
};

/**
 * @brief A class for dual structure \f$\mathrm B(e,e,n+1)\f$ canonical factors.
 *
 * Their structure is described in Bessis, Corran, _Non-Crossing Partitions of
 * Type_ \f$(e, e, r)\f$, 2004, arXiv: [math/0101158
 * [math.GR]](https://arxiv.org/abs/math/0403400).
 *
 * They are most of the time thought of as matrices (elements of semi-direct
 * product \f$\mathrm G(e,e,n+1)=\Delta(e,e,n+1)\mathop{⋊}\mathfrak S_{n+1}\f$,
 * where \f$\Delta(e,e,n+1)\f$ is the group of diagonal matrices of dimension
 * \f$n+1\f$ whose determinant is \f$1\f$ and whose diagonal coefficients are in
 * \f$\mathbb U_e\f$, and where \f$\mathfrak S_{n+1}\f$ is seen as the group of
 * permutation matrices).
 *
 * In practice, only \f$\mathrm O(n)\f$ memory is needed as these matrices are very
 * sparse. They are represented using two vectors, one for the permutation and
 * one for the coefficients (_i.e._ the diagonal part of the semi-direct
 * product).
 *
 * Let \f$\zeta_{ne}=\mathrm e^{i\frac{1}{en}\tau}\f$ and \f$\zeta_{e}=\mathrm
 * e^{i\frac{1}{e}\tau}\f$. To an element in \f$\mathrm G(e,e,n+1)\f$ one can
 * associate a non-crossing partition on \f$\mathbb U_{ne}\cup\{0\}\f$ (see the
 * article). To do this, it is more convenient to consider elements as
 * coefficiented permutations of \f$0, \zeta_{ne}^{-1}, \zeta_{ne}^{-2},\ldots,
 * \zeta_{ne}^{-n}\f$. Therefore, we will often identify \f$k\in[\![1,n]\!]\f$
 * with \f$\zeta_{ne}^{-k}\f$.
 */
class Underlying {

  public:
    /**
     * @brief Parameter type.
     */
    using Parameter = EENParameter;

  private:
    /**
     * @brief The group parameter.
     *
     * An instance with members `e` and `n` set to \f$e\f$ and \f$n\f$
     * respectively is the parameter of \f$\mathrm B(e,e,n+1)\f$.
     */
    Parameter een_index;

    /**
     * @brief The induced permutation.
     *
     * The inverse table of the permutation induced by the matrix on indexes
     * (from \f$0\f$ to \f$n\f$), through its action by left multiplication
     * on vectors.
     *
     * This permutation sends \f$i\f$ to the \f$j\f$ such that the \f$i\f$-th
     * coordinate is the \f$j\f$-th after multiplication. That is to say,
     * \f$j\f$ is the index such that the coefficient at \f$(j, i)\f$ is non
     * zero. Thus, in the inverse table, the integer \f$k\f$ at the
     * \f$f\f$-position is the one such that the coefficient at \f$(i,k)\f$ is
     * non zero.
     */
    std::vector<i16> permutation_table;

    /**
     * @brief The multiplicating coefficients.
     *
     * The multiplicating coefficients. These are \f$e\f$-th roots of unity,
     * with the added condition of their product's being one. They are
     * represented by integer ranging between \f$0\f$ and \f$e - 1\f$, with
     * \f$i\f$ standing for \f$\zeta_e^i\f$.
     *
     * `coefficient_table` is the diagonal coefficient matrix, when writing the
     * matrix as a profuct of a diagonal matrix on the left and a permutation
     * matrix on the right.
     */
    std::vector<i16> coefficient_table;

  public:
    /**
     * @brief Maximum value for the \f$n\f$ index.
     *
     * The greatest \f$n\f$ value that may be used for factors.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const i16 MAX_N_PARAMETER = 256;

    /**
     * @brief Maximum value for the \f$e\f$ index.
     *
     * Not exactly. This is not a strict bound on \f$e\f$. Rather, the
     * real condition is \f$en\f$ should be at most MAX_E_PARAMETER *
     * MAX_N_PARAMETER`.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const i16 MAX_E_PARAMETER = 1;

    /**
     * @brief Converts a string to a parameter.
     *
     * Converts a string to a parameter by tring to parse it as a couple of
     * integers.
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
     * That is to say, the \f$(e,n)\f$ in \f$\mathrm B(e,e,n+1)\f$.
     *
     * @return The parameter.
     */
    inline Parameter get_parameter() const { return een_index; }

    /**
     * @brief Height of the lattice.
     *
     * (_I.e._ the length of
     * the Garside element \f$\Delta\f$ as a word in the generators, which is
     * \f$n+1\f$.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const { return get_parameter().n + 1; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `p` as its
     * parameter.
     *
     * Its tables will have length `p.n + 1`, and will be filled with
     * zeros (thus this is not a valid factor). It should be initialized it with
     * `identity()`, `delta()`, or another similar member.
     *
     * @param p The parameter of the factor (also the length of
     * its tables, minus one).
     */
    Underlying(Parameter p);

    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Prints internal representation in `os`.
     *
     * Prints private members `een_index`, `permutation_table` and
     * `coefficient_table`, typically for debugging.
     *
     * @param os The output stream it is printed in.
     */
    void debug(IndentedOStream &os) const;

    /**
     * @brief Computes the partition associated with the factor.
     *
     * It is then written in `i16` array `x`.
     *
     * A partition is represented by an integer array,
     * with the following convention: given an array \f$\mathrm A\f$
     * representing a partition, then for all \f$i\in[\![0,en]\!]\f$, \f$\mathrm
     * A[i]\f$ is the minimum of the cell \f$i\f$ belongs to.
     *
     * Runs in \f$\mathrm O(en)\f$ time.
     *
     * See Bessis, Corran, _Non-Crossing Partitions of
     * Type_ \f$(e, e, r)\f$, 2004, arXiv: [math/0101158
     * [math.GR]](https://arxiv.org/abs/math/0403400).
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
     * representing a partition, then for all \f$i\in[\![0,en]\!]\f$, \f$\mathrm
     * A[i]\f$ is the minimum of the cell \f$i\f$ belongs to.
     *
     * Runs in \f$\mathrm O(en)\f$ time.
     *
     * See Bessis, Corran, _Non-Crossing Partitions of
     * Type_ \f$(e, e, r)\f$, 2004, arXiv: [math/0101158
     * [math.GR]](https://arxiv.org/abs/math/0403400).
     *
     * @param x A `i16` array.
     */
    void of_partition(const i16 *x);

    /**
     * @brief Prints the factor to `os`.
     *
     * It is printed as a product in the Bessis-Corran short and long
     * generators, in a way that spells the disjunct cycle decomposition.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    /**
     * @brief Sets the factor to the identity.
     *
     * (_I.e._ sets the permutation table to the one of the identity
     * matrix.)
     *
     * Linear in \f$n\f$.
     */
    void identity();

    /**
     * @brief Sets the factor to the Garside element.
     *
     * (_I.e._ sets the tables to the ones of \f$\Delta\f$. \f$\Delta\f$ is the
     * matrix such that \f$[\Delta]_{0,0}=\zeta_e^{-1}\f$,
     * \f$[\Delta]_{n,1}=\zeta_e\f$, and then for \f$i\in[\![1,n-1]\!]\f$,
     * \f$[\Delta]_{i,i+1}=1\f$.)
     *
     * Linear in \f$n\f$.
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
     * Runs in \f$\mathrm O(en)\f$ time, although it uses a `thread_local`
     * square matrix of dimension `MAX_N_PARAMETER * MAX_E_PARAMETER + 1`.
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
     * Runs in \f$\mathrm O(en)\f$ time, although it uses a `thread_local`
     * square matrix of dimension `MAX_N_PARAMETER * MAX_E_PARAMETER + 1`.
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
     * they represent the same matrix).
     *
     * Linear in \f$n\f$.
     *
     * @param b Second operand.
     * @return If `*this` and `b` are equal.
     */
    inline bool compare(const Underlying &b) const {
        return (permutation_table == b.permutation_table) &&
               (coefficient_table == b.coefficient_table);
    }

    /**
     * @brief Product computations.
     *
     * Computes the product of `*this` and `b` (_i.e._ matrix product), under
     * the assumption that it is a factor.
     *
     * Linear in \f$n\f$.
     *
     * @param b Second (right) operand.
     * @return The product of `*this` and `b`.
     */
    Underlying product(const Underlying &b) const;

    /**
     * @brief Left complement computations.
     *
     * Computes the left complement of `*this` to `b`, under the assumption
     * that `*this` divides `b`.
     *
     * Linear in \f$n\f$.
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
     * Linear in \f$n\f$.
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
    void randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (_i.e._ Bessis-Corran short symmetric and
     * assymmetric generators).
     *
     * See Bessis, Corran, _Non-Crossing Partitions of
     * Type_ \f$(e, e, r)\f$, 2004, arXiv: [math/0101158
     * [math.GR]](https://arxiv.org/abs/math/0403400).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by \f$\Delta^k\f$.
     *
     * Linear in \f$n\f$ (in particular, does not depend on `k`).
     *
     * @param k The exponent.
     */
    void delta_conjugate_mut(i16 k);

    /**
     * @brief Hashes the factor.
     *
     * Linear in \f$n\f$.
     *
     * @return The hash.
     */
    std::size_t hash() const;

  private:
    /**
     * @brief Computes the factor associated to the inverse of the matrix
     * of this factor.
     *
     * Linear in \f$n\f$.
     *
     * @return The factor associated to the inverse of the matrix of
     * `*this`.
     */
    Underlying inverse() const;
};

/**
 * @brief Class for dual Garside structure \f$\mathrm B(e,e,n+1)\f$ complex
 * braid groups canonical factors.
 */
using Factor = FactorTemplate<Underlying>;

/**
 * @brief Class for dual Garside structure \f$\mathrm B(e,e,n+1)\f$ complex
 * braid groups elements.
 */
using Braid = BraidTemplate<Factor>;

} // namespace dual_complex

/**
 * @brief Inserts a parameter in the output stream.
 *
 * Syntactic sugar for `dual_complex::EENParameter.print()`.
 *
 * @param p The parameeter to be inserted.
 * @return A reference to `*this`, so that `<<` may be chained.
 */
template <>
inline IndentedOStream &
    IndentedOStream::operator<< <dual_complex::EENParameter>(
        const dual_complex::EENParameter &p) {
    p.print(*this);
    return *this;
}

} // namespace garcide
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
     * @brief `e` field.
     *
     * It is the `e` such that coefficients are in \f$\mathbb U_e\f$.
     */
    i16 e;

    /**
     * @brief `n` field.
     *
     * Factors are (conceptually) square matrix of dimension `n+1`.
     */
    i16 n;

    /**
     * @brief Construct a new `EENParameter`.
     *
     * @param e Its `e` field.
     * @param n Its `n` field.
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
 * Type \f$(e, e, r)\f$_, 2004, arXiv: [math/0101158
 * [math.GR]](https://arxiv.org/abs/math/0403400).
 *
 * They are most of the time thought of as matrix (elements of semi-direct
 * product \f$\Delta(e,e,n+1)\mathop{⋊}\mathfrak S_{n+1}\f$, where
 * \f$\Delta(e,e,n+1)\f$ is the group of diagonal matrices of dimension
 * \f$n+1\f$ whose determinant is \f$1\f$ and whose diagonal coefficients are in
 * \f$\mathbb U_e\f$, and where \f$\mathfrak S_{n+1}\f$ is seen as the group of
 * permutation matrices).
 */
class Underlying {

  public:
    /**
     * @brief Parameter type.
     */
    using Parameter = EENParameter;

  private:
    Parameter een_index;

    // The induced permutation, where 0 is the 0-th coordinates, and then the
    // next n coordinates represent the powers w ^ -i (with i ranging between 0
    // and n - 1) of an ne-th root of unity w. We use the same conventions as
    // before: letting sigma be the induced permutation, then
    // permutation_table[i] is sigma^(-1)(i).
    std::vector<i16> permutation_table;

    // The multiplicating coefficients. These are e-th roots of unity, with the
    // added condition of their product's being one. They are represented by
    // integer ranging between 0 and e - 1, with i standing for w^(ei) (with w
    // the same root as for the permutation_table).
    std::vector<i16> coefficient_table;

  public:
    /**
     * @brief Maximum value for `een_index.n`.
     *
     * The greatest `een_index.n` value that may be used for braids.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const i16 MAX_N_PARAMETER = 256;

    /**
     * @brief Maximum value for `een_index.e`.
     *
     * Not exactly. This is not a strict bound on `een_index.e`. The
     * real condition is `een_index.e * een_index.n <=
     * MAX_E_PARAMETER * MAX_N_PARAMETER`.
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
     * (_I.e._, letting `n` be the number of of strands, the length of
     * the Garside element \f$\Delta\f$ as a word in the generators, which is
     * \f$n-1\f$.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const { return get_parameter().n + 1; }

    // Constructor
    Underlying(Parameter p);

    void of_string(const std::string &str, size_t &pos);

    void debug(IndentedOStream &os) const;

    void assign_partition(i16 *x) const;

    void of_partition(const i16 *x);

    // print to os. Be wary, as it side-effects!
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    Underlying left_meet(const Underlying &b) const;

    inline Underlying right_meet(const Underlying &b) const {
        return left_meet(b);
    }

    // Equality check.
    // We check whether the underlying permutation table are (pointwise) equal.
    bool compare(const Underlying &b) const;

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
     * Type \f$(e, e, r)\f$_, 2004, arXiv: [math/0101158
     * [math.GR]](https://arxiv.org/abs/math/0403400).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    // Conjugate by delta^k.
    // Used to speed up calculations compared to the default implementation.
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
     * @brief Computes the factor associated with the inverse of the matrix
     * of this factor.
     *
     * Linear in the parameter.
     *
     * @return The factor associated with the inverse of the matrix of
     * `*this`.
     */
    Underlying inverse() const;
};

using Factor = FactorTemplate<Underlying>;

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
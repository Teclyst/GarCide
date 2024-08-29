/**
 * @file standard_complex.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for \f$\mathrm B(e, e, n)\f$ groups (semi-classic Garside
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
 * @brief Namespace for complex braid groups \f$\mathrm B(e,e,n)\f$,
 * semi-classic Garside structure.
 */
namespace standard_complex {

/**
 * @brief Class for complex braid groups parameters.
 *
 * An instance with members `e` and `n` set to \f$e\f$ and \f$n\f$ respectively
 * is the parameter of \f$\mathrm B(e,e,n)\f$.
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
     * Factors are (conceptually) square matrix of dimension `n`.
     */
    i16 n;

    /**
     * @brief Construct a new `EENParameter`.
     *
     * @param e Its `e` member.
     * @param n Its `n` member.
     */
    EENParameter(i16 e, i16 n) : e(e), n(n) {}

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
 * @brief A class for \f$\mathrm B(e,e,n)\f$ semi-classic canonical factors.
 *
 * Their structure is described in Corran, Picantin, _A New Garside Structure
 * for the Braid Groups of Type_ \f$(e, e, r)\f$, 2011, arXiv: [arxiv/0901.0645
 * [math.GR]](https://arxiv.org/abs/0901.0645) and Neaime, _Interval Garside
 * Structures for the Complex Braid Groups_ \f$B(e,e,n)\f$, 2019, arXiv:
 * [arXiv:1707.06864 [math.GR]](https://arxiv.org/abs/1707.06864).
 *
 * They are most of the time thought of as matrices (elements of semi-direct
 * product \f$\mathrm G(e,e,n)=\Delta(e,e,n)\mathop{⋊}\mathfrak S_{n}\f$,
 * where \f$\Delta(e,e,n)\f$ is the group of diagonal matrices of dimension
 * \f$n\f$ whose determinant is \f$1\f$ and whose diagonal coefficients are in
 * \f$\mathbb U_e\f$, and where \f$\mathfrak S_{n}\f$ is seen as the group of
 * permutation matrices).
 *
 * In practice, only \f$\mathrm O(n)\f$ memory is needed as these matrices are very
 * sparse. They are represented using two vectors, one for the permutation and
 * one for the coefficients (_i.e._ the diagonal part of the semi-direct
 * product).
 *
 * Let \f$\zeta_{e}=\mathrm e^{i\frac{1}{e}\tau}\f$.
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
     * respectively is the parameter of \f$\mathrm B(e,e,n)\f$.
     */
    Parameter een_index;

    /**
     * @brief The induced permutation.
     *
     * The inverse table of the permutation induced by the matrix on indexes
     * (from \f$0\f$ to \f$n - 1\f$), through its action by left multiplication
     * on vectors.
     *
     * This permutation sends \f$i\f$ to the \f$j\f$ such that the \f$i\f$-th
     * coordinate is the \f$j\f$-th after multiplication. That is to say,
     * \f$j\f$ is the index such that the coefficient at \f$(j, i)\f$ is non
     * zero. Thus, in the inverse table, the integer \f$k\f$ at the
     * \f$f\f$-position is the one such that the coefficient at \f$(i,k)\f$ is
     * non zero (which is \f$c_i\f$ in George Neaime, _Interval Garside
     * Structures for the Complex Braid Groups_ \f$B(e,e,n)\f$,
     * [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
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
     * That is to say, the \f$(e,n)\f$ in \f$\mathrm B(e,e,n)\f$.
     *
     * @return The parameter.
     */
    inline Parameter get_parameter() const { return een_index; }

    /**
     * @brief Height of the lattice.
     *
     * (_I.e._, letting `n` be the number of of strands, the length of
     * the Garside element \f$\Delta\f$ as a word in the generators, which is
     * \f$n(n+1)\f$.)
     *
     * @return The height of the lattice.
     */
    inline i16 lattice_height() const {
        return get_parameter().n * (get_parameter().n + 1);
    }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `p` as its
     * parameter.
     *
     * Its tables will have length `p.n`, and will be filled with
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
     * @brief Prints the factor to `os`.
     *
     * Prints to `os` Neaime's shortest word representative in the
     * Corran-Picantin generators (See George Neaime, _Interval Garside
     * Structures for the Complex Braid Groups_ \f$B(e,e,n)\f$,
     * [arXiv:1707.06864 [math.GR]](https://arxiv.org/abs/1707.06864))
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
     * (_I.e._ sets the tables to the ones of \f$\lambda_{en}\f$.
     * \f$\lambda_{en}\f$ is the matrix such that
     * \f$[\lambda_{en}]_{0,0}=\zeta_e^{-n+1}\f$, and then
     * for \f$i\in[\![1,n-1]\!]\f$, \f$[\Delta]_{i,i}=\zeta_e\f$.)
     *
     * Linear in  \f$n\f$.
     */
    void delta();

    /**
     * @brief Computes the left meet of `*this` and `b`.
     *
     * This is done by extracting it atom by atom in Neaime's normal form for
     * factors.
     *
     * Runs in \f$\mathrm O(n^2+en)\f$ time.
     *
     * @param b Second operand.
     * @return The meet of `*this` and `b`.
     */
    Underlying left_meet(const Underlying &b) const;

    /**
     * @brief Computes the right meet of `*this` and `b`.
     *
     * Uses an identity to reduce the calculation to the left case.
     *
     * @param b Second operand.
     * @return The right meet of `*this` and `b`.
     */
    inline Underlying right_meet(const Underlying &b) const {
        return inverse().left_meet(b.inverse()).inverse();
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
     * Returns the list of the atoms (_i.e._ Corran-Picantin \f$s\f$ and \f$t\f$
     * generators).
     *
     * See Corran, Picantin, _A New Garside Structure
     * for the Braid Groups of Type_ \f$(e, e, r)\f$, 2011, arXiv:
     * [arxiv/0901.0645 [math.GR]](https://arxiv.org/abs/0901.0645)<.
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
     * @brief Gets the direct permutation table associated with the factor.
     *
     * (I.e., perm[i] is the image of i by the permutation.)
     * Used as dividing by atoms requires columns to be easily findable.
     *
     * @param dir_perm The table that is to be filled.
     */
    void direct(i16 *dir_perm) const;

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

    // Checks if s_i left divides `*this`. (See George Neaime, Interval Garside
    // Structures for the Complex Braid Groups, Proposition 3.13,
    // [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor.
    // @param i The integer i for which we check if s_i left divides the factor.
    inline bool is_s_left_divisor(i16 i) const {
        return (permutation_table[i - 1] > permutation_table[i - 2])
                   ? (coefficient_table[i - 1] != 0)
                   : (coefficient_table[i - 2] == 0);
    }

    // Checks if t_i left divides `*this`. (See George Neaime, Interval Garside
    // Structures for the Complex Braid Groups, Proposition 3.13,
    // [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor.
    // @param i The integer i for which we check if t_i left divides the factor.
    inline bool is_t_left_divisor(i16 i) const {
        return (permutation_table[1] > permutation_table[0])
                   ? (coefficient_table[1] != 0)
                   : (coefficient_table[0] ==
                      ((i == 0) ? 0 : get_parameter().e - i));
    }

    // Left multiplies by s_i (or divides, which is the same as it has order 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by s_i the factor.
    inline void s_left_multiply(i16 *dir_perm, i16 i) {
        std::swap(coefficient_table[i - 1], coefficient_table[i - 2]);
        std::swap(permutation_table[i - 1], permutation_table[i - 2]);
        std::swap(dir_perm[permutation_table[i - 1]],
                  dir_perm[permutation_table[i - 2]]);
    }

    // Left multiplies by t_i (or divides, which is the same as it has order 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by t_i the factor.
    inline void t_left_multiply(i16 *dir_perm, i16 i) {
        std::swap(coefficient_table[0], coefficient_table[1]);
        std::swap(permutation_table[0], permutation_table[1]);
        coefficient_table[0] = Rem(coefficient_table[0] - i, get_parameter().e);
        coefficient_table[1] = Rem(coefficient_table[1] + i, get_parameter().e);
        std::swap(dir_perm[permutation_table[0]],
                  dir_perm[permutation_table[1]]);
    }

    // Right multiplies by s_i (or divides, which is the same as it has order
    // 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by s_i the factor.
    inline void s_right_multiply(i16 *dir_perm, i16 i) {
        std::swap(permutation_table[dir_perm[i - 1]],
                  permutation_table[dir_perm[i - 2]]);
        std::swap(dir_perm[i - 1], dir_perm[i - 2]);
    }

    /** Right multiplies by t_i.
     *
     * (Or divides, which is the same as it has order \f$2\f$.)
     *
     * @param dir_perm A table that holds the direct table of the permutation
     * induced by the factor. Modified by this function so that it remains up to
     * date.
     * @param i The integer \f$i\f$ for which we multiply the factor by
     * \f$t_i\f$.
     */
    inline void t_right_multiply(i16 *dir_perm, i16 i) {
        std::swap(permutation_table[dir_perm[0]],
                  permutation_table[dir_perm[1]]);
        coefficient_table[dir_perm[0]] =
            Rem(coefficient_table[dir_perm[0]] - i, get_parameter().e);
        coefficient_table[dir_perm[1]] =
            Rem(coefficient_table[dir_perm[1]] + i, get_parameter().e);
        std::swap(dir_perm[0], dir_perm[1]);
    }
};

/**
 * @brief Class for semi-classic Garside structure \f$\mathrm B(e,e,n+1)\f$ complex
 * braid groups canonical factors.
 */
using Factor = FactorTemplate<Underlying>;

/**
 * @brief Class for semi-classic Garside structure \f$\mathrm B(e,e,n+1)\f$ complex
 * braid groups elements.
 */
using Braid = BraidTemplate<Factor>;

} // namespace standard_complex

/**
 * @brief Inserts a parameter in the output stream.
 *
 * Syntactic sugar for `standard_complex::EENParameter.print()`.
 *
 * @param p The parameter to be inserted.
 * @return A reference to `*this`, so that `<<` may be chained.
 */
template <>
inline IndentedOStream &
    IndentedOStream::operator<< <standard_complex::EENParameter>(
        const standard_complex::EENParameter &p) {
    p.print(*this);
    return *this;
}

} // namespace garcide
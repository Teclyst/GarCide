/**
 * @file standard_complex.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for B(e, e, n) groups (semi-classic Garside structure).
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

#include "garcide.h"

namespace garcide {

namespace standard_complex {

// We represent B(e, e, n).
struct EENParameter {
    sint16 e;
    sint16 n;

    EENParameter(sint16 e, sint16 n) : e(e), n(n) {}

    inline bool compare(const EENParameter &p) const {
        return ((e == p.e) && (n == p.n));
    }

    inline bool operator==(const EENParameter &p) const { return compare(p); }

    inline bool operator!=(const EENParameter &p) const { return !compare(p); }

    void print(IndentedOStream &os) const;
};

class Underlying {

  public:
    using Parameter = EENParameter;

  private:
    Parameter een_index;

    /**
     * @brief The induced permutation.
     *
     * The table of the permutation induced by the matrix on indexes (from 0
     * to `een_index.n` - 1), through its action by left
     * multiplication on vectors.
     *
     * This permutation sends i to the j such that the i-th coordinate is the
     * j-th after multiplication. That is to say, j is the index such that the
     * coefficient at (i, j) is non zero (c_i in George Neaime, Interval Garside
     * Structures for the Complex Braid Groups,
     * [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
     */
    std::vector<sint16> permutation_table;

    /**
     * @brief The multiplicating coefficients.
     *
     * Let e = `een_index.e` and n = `een_index.n`.
     *
     * The multiplicating coefficients. These are e-th roots of unity, with the
     * added condition of their product's being one. They are represented by
     * integer ranging between 0 and e - 1, with i standing for zeta^i (with
     * zeta an e-th root of unity).
     *
     * `coefficient_table` is the diagonal coefficient matrix, when considering
     * G(e, e, n) as the semi-direct product delta(e, e, n) ⋊ S(n).
     */
    std::vector<sint16> coefficient_table;

  public:
    /**
     * @brief Maximum braid index.
     *
     * The greatest index that may be used for braids.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const sint16 MAX_N = 256;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    sint16 lattice_height() const;

    // Constructor
    Underlying(Parameter p);

    void of_string(const std::string &str, size_t &pos);

    // Prints to `os` the internal representation of the factor.
    // @param os The `std::ostream` we are printing to.
    void debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints to `os` Neaime's shortest word representative in the
     * Corran-Picantin generators (See George Neaime, Interval Garside
     * Structures for the Complex Braid Groups,
     * [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)))
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    // Gets the direct permutation table associated with the factor.
    // (i.e., perm[i] is the image of i by the permutation).
    // Used as dividing by atoms requires columns to be easily findable.
    // @param dir_perm The table that is to be filled.
    void direct(sint16 *dir_perm) const;

    // Checks if s_i left divides `*this`. (See George Neaime, Interval Garside
    // Structures for the Complex Braid Groups, Proposition 3.13,
    // [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor.
    // @param i The integer i for which we check if s_i left divides the factor.
    inline bool is_s_left_divisor(const sint16 *dir_perm, sint16 i) const {
        return (permutation_table[i - 1] > permutation_table[i - 2])
                   ? (coefficient_table[i - 1] != 0)
                   : (coefficient_table[i - 2] == 0);
    };

    // Checks if t_i left divides `*this`. (See George Neaime, Interval Garside
    // Structures for the Complex Braid Groups, Proposition 3.13,
    // [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor.
    // @param i The integer i for which we check if t_i left divides the factor.
    inline bool is_t_left_divisor(const sint16 *dir_perm, sint16 i) const {
        return (permutation_table[1] > permutation_table[0])
                   ? (coefficient_table[1] != 0)
                   : (coefficient_table[0] ==
                      ((i == 0) ? 0 : get_parameter().e - i));
    };

    // Left multiplies by s_i (or divides, which is the same as it has order 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by s_i the factor.
    inline void s_left_multiply(sint16 *dir_perm, sint16 i) {
        std::swap(coefficient_table[i - 1], coefficient_table[i - 2]);
        std::swap(permutation_table[i - 1], permutation_table[i - 2]);
        std::swap(dir_perm[permutation_table[i - 1]],
                  dir_perm[permutation_table[i - 2]]);
    };

    // Left multiplies by t_i (or divides, which is the same as it has order 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by t_i the factor.
    inline void t_left_multiply(sint16 *dir_perm, sint16 i) {
        std::swap(coefficient_table[0], coefficient_table[1]);
        std::swap(permutation_table[0], permutation_table[1]);
        coefficient_table[0] = Rem(coefficient_table[0] - i, get_parameter().e);
        coefficient_table[1] = Rem(coefficient_table[1] + i, get_parameter().e);
        std::swap(dir_perm[permutation_table[0]],
                  dir_perm[permutation_table[1]]);
    };

    // Right multiplies by s_i (or divides, which is the same as it has order
    // 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by s_i the factor.
    inline void s_right_multiply(sint16 *dir_perm, sint16 i) {
        std::swap(permutation_table[dir_perm[i - 1]],
                  permutation_table[dir_perm[i - 2]]);
        std::swap(dir_perm[i - 1], dir_perm[i - 2]);
    };

    // Right multiplies by t_i (or divides, which is the same as it has order
    // 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by t_i the factor.
    inline void t_right_multiply(sint16 *dir_perm, sint16 i) {
        std::swap(permutation_table[dir_perm[0]],
                  permutation_table[dir_perm[1]]);
        coefficient_table[dir_perm[0]] =
            Rem(coefficient_table[dir_perm[0]] - i, get_parameter().e);
        coefficient_table[dir_perm[1]] =
            Rem(coefficient_table[dir_perm[1]] + i, get_parameter().e);
        std::swap(dir_perm[0], dir_perm[1]);
    };

    Underlying left_meet(const Underlying &b) const;

    inline Underlying right_meet(const Underlying &b) const {
        return inverse().left_meet(b.inverse()).inverse();
    };

    // Equality check.
    // We check whether the underlying permutation table are (pointwise) equal.
    bool compare(const Underlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    Underlying inverse() const;

    // product under the hypothesis that it is still simple.
    Underlying product(const Underlying &b) const;

    // Under the assumption a <= b, a.left_complement(b) computes
    // The factor c such that ac = b.
    Underlying left_complement(const Underlying &b) const;

    Underlying right_complement(const Underlying &b) const;

    // Generate a random factor.
    void randomize();

    // List of atoms.
    std::vector<Underlying> atoms() const;

    // Conjugate by delta^k.
    // Used to speed up calculations compared to the default implementation.
    void delta_conjugate_mut(sint16 k);

    std::size_t hash() const {
        std::size_t h = 0;
        for (sint16 i = 0; i < get_parameter().n; i++) {
            h = h * 31 + permutation_table[i];
        }
        for (sint16 i = 0; i < get_parameter().n; i++) {
            h = h * 31 + coefficient_table[i];
        }
        return h;
    }
};

using Factor = FactorTemplate<Underlying>;

using Braid = BraidTemplate<Factor>;

} // namespace standard_complex

template <>
IndentedOStream &IndentedOStream::operator<< <standard_complex::EENParameter>(
    const standard_complex::EENParameter &p);

} // namespace cgarside
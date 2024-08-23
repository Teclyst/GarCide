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

/**
 * @brief A class for \f$\mathrm B(e,e,n+1)\f$ dual canonical factors.
 *
 * Their structure is described in Bessis, Corran, _Non-crossing partitions of
 * type \f$(e, e, r)\f$_, 2004, arXiv: [math/0101158
 * [math.GR]](https://arxiv.org/abs/math/0403400).
 */
class Underlying {

  public:
    using Parameter = EENParameter;

  private:
    Parameter een_index;

    // The induced permutation, where 0 is the 0-th coordinates, and then the
    // next n coordinates represent the powers w ^ -i (with i ranging between 0
    // and n - 1) of an ne-th root of unity w. We use the same conventions as
    // before: letting sigma be the induced permutation, then
    // permutation_table[i] is sigma^(-1)(i).
    std::vector<sint16> permutation_table;

    // The multiplicating coefficients. These are e-th roots of unity, with the
    // added condition of their product's being one. They are represented by
    // integer ranging between 0 and e - 1, with i standing for w^(ei) (with w
    // the same root as for the permutation_table).
    std::vector<sint16> coefficient_table;

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
    static const sint16 MAX_N = 256;

    /**
     * @brief Maximum value for `een_index.e`.
     *
     * Not exactly. This is not a strict bound on `een_index.e`. The
     * real condition is `een_index.e * een_index.n <=
     * MAX_E * MAX_N`.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const sint16 MAX_E = 1;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    sint16 lattice_height() const;

    // Constructor
    Underlying(Parameter p);

    void of_string(const std::string &str, size_t &pos);

    void debug(IndentedOStream &os) const;

    void assign_partition(sint16 *x) const;

    void of_partition(const sint16 *x);

    // print to os. Be wary, as it side-effects!
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    Underlying left_meet(const Underlying &b) const;

    inline Underlying right_meet(const Underlying &b) const {
        return left_meet(b);
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

    std::size_t hash() const;
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

} // namespace dual_complex

template <>
IndentedOStream &IndentedOStream::operator<< <dual_complex::EENParameter>(
    const dual_complex::EENParameter &p);

} // namespace garcide
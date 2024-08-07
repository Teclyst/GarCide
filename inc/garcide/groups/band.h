/**
 * @file band.h
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

namespace garcide::band {

class Underlying {
  public:
    using Parameter = sint16;

  private:
    Parameter number_of_strands;

    std::vector<sint16> permutation_table;

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
    static const Parameter MAX_NUMBER_OF_STRANDS = 256;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    sint16 lattice_height() const;

    // Constructor
    Underlying(sint16 n);

    sint16 at(size_t i) const { return permutation_table[i]; }

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting `W = (\s | \t)*` be the language of whitespaces and `Z = -? ([1 -
     * 9] [0 - 9]* | 0)` be the language of integers, accepted strings are those
     * represented by the regular expression `\(W Z W,? W Z W\)`, under the
     * additional hypothesis that the two integers lie in [`1`, `Parameter`] and
     * are not equal.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches `\(W Z W,? W Z W\)`, or if there is one, if
     * either integer does not belong to [`1`, `Parameter`], or both are equal.
     */
    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Prints internal representation to `os`.
     *
     * Prints the factor's `permutation_table` to `os`, typically for debugging
     * purposes.
     *
     * @param os The output stream it prints to.
     */
    void debug(IndentedOStream &os) const;

    void assign_partition(sint16 *x) const;

    void of_partition(const sint16 *x);

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints the factor to `os` as a product of atoms.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    Underlying left_meet(const Underlying &b) const;

    Underlying right_meet(const Underlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
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

    size_t hash() const;

    void of_ballot_sequence(const sint8 *s);
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

#ifdef USE_CLN

void ballot_sequence(sint16 n, cln::cl_I k, sint8 *s);

const cln::cl_I &get_catalan_number(sint16 n);

#endif

} // namespace garcide::band

#endif
/**
 * @file octahedral.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for \f$\mathbf B\f$-series Artin groups (dual Garside structure).
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
 * @brief Namespace for \f$\mathbf B\f$-series Artin groups, dual Garside structure.
 */
namespace garcide::octahedral {

class Underlying {

  protected:
    i16 PresentationParameter;

    std::vector<i16> permutation_table;

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
    static const i16 MaxBraidIndex = 256;

    typedef i16 Parameter;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    i16 lattice_height() const;

    // Constructor
    Underlying(i16 n);

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting `W = (\s | \t)*` be the language of whitespaces and `Z = -? ([1 -
     * 9] [0 - 9]* | 0)` be the language of integers, accepted strings are those
     * represented by regular expression `\(W Z W,? W Z W\) | Z`, under the
     * additional hypothesis that in the first case the two integers are not
     * equal mod `PresentationParameter'.
     *
     * The first case stands for short generators (double, symmetric,
     * transpositions), and the second one for long generators (transposition of
     * antipodals points).
     * 
     * "D" is also accepted and is parsed as \f$\Delta\f$.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception `InvalidStringError`: Thrown when there is no subword starting
     * from `pos` that matches `\(W Z W,? W Z W\) | Z | D`, or if there is one, it
     * matches `\(W Z W,? W Z W\)`, and both integers are equal mod
     * `PresentationParameter'.
     */
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
    void delta_conjugate_mut(i16 k);

    std::size_t hash() const {
        std::size_t h = 0;
        for (i16 i = 1; i <= 2 * get_parameter(); i++) {
            h = h * 31 + permutation_table[i];
        }
        return h;
    }
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

} // namespace CGarside

#endif
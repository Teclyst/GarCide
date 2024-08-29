/**
 * @file dihedral.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file for euclidian lattices \f$\mathbf I\f$-series Artin groups
 * (dual Garside structure).
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

#include "garcide/garcide.hpp"

/**
 * @brief Namespace for \f$\mathbf I\f$-series Artin groups, dual Garside
 * structure.
 */
namespace garcide::dihedral {

/**
 * @brief Exception thrown when an illegal operation is attempted.
 *
 * In practice, complements where a factor that should divide the other does
 * not, or products whose result is not simple.
 */
struct NotBelow {};

/**
 * @brief A class for dual Garside structure \f$\mathbf I\f$-series Artin groups
 * canonical factors.
 *
 * They are represented by elements of dihedral group \f$\mathrm D_{2n}\f$,
 * using constant space (as a function of \f$n\f$).
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
     * @brief Group parameter.
     *
     * It is the \f$n\f$ in \f$\mathrm D_{2n}\f$: the number of vertices in the
     * regular polygon \f$\mathrm D_{2n}\f$ naturally acts.
     */
    Parameter number_of_vertices;

    /**
     * @brief Factor type.
     *
     * `0` for identity, `1` for a rotation, `2` for a reflection.
     *
     * In the dual structure, only one rotation is a factor, thus `1`
     * also stands for the Garside element \f$\Delta\f$.
     */
    i16 type;

    /**
     * @brief The vertex the factor sends \f$0\f$ on.
     *
     * Only relevant if `type` is `2` (_i.e_ for reflections). Otherwise, the
     * type already determines the factor.
     */
    i16 vertex;

  public:
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
     * That is to say, the \f$n\f$ in \f$\mathrm D_{2n}\f$.
     *
     * @return The group parameter.
     */
    inline Parameter get_parameter() const { return number_of_vertices; }

    /** @brief Height of the lattice.
     *
     * (_I.e._, letting `n` be the parameter, the length of
     * \f$\Delta\f$ as a word in the generators, which is
     * \f$2\f$.)
     *
     * @return The height of the lattice.
     */
    i16 lattice_height() const { return 2; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `n` as its
     * parameter.
     *
     * It will be initialized as the identity factor.
     *
     * @param n The parameter of the factor.
     */
    Underlying(i16 n);

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting \f$Z = \texttt{-}? ([\texttt{1} - \texttt{9}] [\texttt{0} -
     * \texttt{9}]^* \mid \texttt{0})\f$ be a regular expression representing
     * integers, accepted strings are those represented by regular expression
     * \f$(\texttt{s} \texttt{_}?)? Z\mid \texttt{D}\f$, ignoring whitespaces.
     * The integer is taken modulo the parameter \f$n\f$.
     *
     * Considering a circular ordering on the vertices of the regular
     * \f$n\f$-gon, \f$\texttt{s_}i\f$, \f$\texttt{s}i\f$ and \f$i\f$ all stand
     * for the reflection \f$\s_i\f$ that sends vertex \f$0\f$ to vertex \f$i\f$
     * (once again, mod \f$n\f$).
     *
     * \f$\texttt{D}\f$ represents \f$\Delta\f$.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception InvalidStringError Thrown when there is no subword starting
     * from `pos` that matches the expression.
     */
    void of_string(const std::string &str, size_t &pos);

    /**
     * @brief Prints internal representation in `os`.
     *
     * Prints private member `number_of_vertices`, `type` and `vertex`,
     * typically for debugging.
     *
     * @param os The output stream it is printed in.
     */
    void debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor to `os`.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    /**
     * @brief Sets the factor to the identity.
     *
     * Constant time.
     */
    inline void identity() { type = 0; }

    /**
     * @brief Sets the factor to the Garside element.
     *
     * Constant time.
     */
    inline void delta() { type = 1; }

    /**
     * @brief Computes the meet of `*this` and `b`.
     *
     * For the dual structure, the left and right meets are equal.
     *
     * Constant time.
     *
     * @param b Second operand.
     * @return The meet of `*this` and `b`.
     */
    Underlying left_meet(const Underlying &b) const;

    /**
     * @brief Computes the meet of `*this` and `b`.
     *
     * For the dual structure, the left and right meets are equal.
     *
     * Constant time.
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
     * they represent the same element in the dihedral group).
     *
     * Linear in the parameter.
     *
     * @param b Second operand.
     * @return If `*this` and `b` are equal.
     */
    inline bool compare(const Underlying &b) const {
        return ((type == b.type) && ((type != 2) || (vertex == b.vertex)));
    }

    /**
     * @brief Product computations.
     *
     * Computes the product of `*this` and `b` (_i.e._ composition of
     * isometries), under the assumption that it is a factor.
     *
     * Constant time.
     *
     * @param b Second (right) operand.
     * @return The product of `*this` and `b`.
     * @exception NotBelow If the product is not simple.
     */
    Underlying product(const Underlying &b) const;

    /**
     * @brief Left complement computations.
     *
     * Computes the left complement of `*this` to `b`, under the assumption that
     * `*this` divides `b`.
     *
     * Constant time.
     *
     * @param b Second operand.
     * @return The left complement of `*this` to `b`.
     * @exception NotBelow If `*this` does not divide `b`.
     */
    Underlying left_complement(const Underlying &b) const;

    /**
     * @brief Right complement computations.
     *
     * Computes the left complement of `*this` to `b`, under the assumption that
     * `*this` divides `b`.
     *
     * Constant time.
     *
     * @param b Second operand.
     * @return The left complement of `*this` to `b`.
     * @exception NotBelow If `*this` does not divide `b`.
     */
    Underlying right_complement(const Underlying &b) const;

    /**
     * @brief Sets `*this` to a random factor.
     */
    void randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (_i.e._ all the reflections).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> atoms() const;

    /**
     * @brief Conjugates by \f$\Delta^k\f$.
     *
     * Constant time.
     *
     * @param k The exponent.
     */
    void delta_conjugate_mut(i16 k);

    /**
     * @brief Hashes the factor.
     * 
     * Constant time.
     * 
     * @return The hash. 
     */
    inline std::size_t hash() const { return (size_t)vertex; }
};

/**
 * @brief Class for dual Garside structure \f$\mathbf I\f$-series Artin groups
 * canonical factors.
 */
using Factor = FactorTemplate<Underlying>;

/**
 * @brief Class for dual Garside structure \f$\mathbf I\f$-series Artin groups elements.
 */
using Braid = BraidTemplate<Factor>;

} // namespace garcide::dihedral
/**
 * @file centralizer.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header (and implementation) file for centralizer computations.
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

#ifndef CENTRALIZER
#define CENTRALIZER

#include "garcide/ultra_summit.h"

/**
 * @brief Namespace for centralizers.
 */
namespace garcide::centralizer {

/**
 * @brief A class representing a collection of generators for a centralizer.
 *
 * A class representing a collection of generators for a centralizer.
 * This is basically a wrapper for a `std::unordered_set`.
 *
 * @tparam B A class representing braids.
 */
template <class B> class Centralizer {
  private:
    /**
     * @brief A set containing generators.
     *
     * The `std::unordered_set` containing generators `*this` is a wrapper of.
     */
    std::unordered_set<B> generators;

  public:
    /**
     * @brief `const` iterator.
     *
     * `const` iterator, that provides read-only iteration on the generators.
     */
    using ConstIterator = typename std::unordered_set<B>::const_iterator;

    /**
     * @brief Iterator to the first element of `*this`.
     *
     * Returns a `const` iterator to the first element of centralizer `*this`.
     *
     * @return An iterator to the first element of `*this`.
     */
    inline ConstIterator begin() const { return generators.begin(); }

    /**
     * @brief Iterator to the after-last element of `*this`.
     *
     * Returns a `const` iterator to the after-last element of centralizer
     * `*this`.
     *
     * @return An iterator to the after-last element of `*this`.
     */
    inline ConstIterator end() const { return generators.end(); }

    /**
     * @brief Insert `b` into `*this`.
     *
     * Insert braid `b` into centralizer `*this`.
     *
     * @param b The braid to be inserted.
     */
    inline void insert(const B &b) { generators.insert(b); }

    /**
     * @brief Returns the number of generators already in `*this`.
     *
     * Returns the number of generators already in centralizer `*this` (_i.e._
     * the size of `*this`).
     *
     * @return The number of generators already in `*this`.
     */
    inline size_t number_of_generators() { return generators.size(); }

    /**
     * @brief Look if `b` is in `*this`.
     *
     * Look if braid `b` is stocked as a generator of centralizer `*this`.
     *
     * @param b The braid looked for.
     * @return if `b` is in `*this`.
     */
    inline bool mem(const B &b) {
        return generators.find(b) != generators.end();
    }

    /**
     * @brief Prints internal representation to `os`.
     *
     * Prints centralizer `*this`'s `generators` to output stream `os`,
     * typically for debugging purposes.
     *
     * @param os The output stream it prints to.
     */
    void debug(IndentedOStream &os = ind_cout) {
        bool is_first = true;
        os << "{   ";
        os.Indent(4);
        os << "generators:";
        os.Indent(4);
        os << EndLine();
        os << "{   ";
        os.Indent(4);
        for (ConstIterator it = begin(); it != end(); it++) {
            if (!is_first) {
                os << "," << EndLine();
            } else {
                is_first = false;
            }
            (*it).debug(os);
        }
        os.Indent(-4);
        os << EndLine();
        os << "}";
        os.Indent(-8);
        os << EndLine();
        os << "}";
    }

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints the centralizer to output stream `os` as a product of atoms.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os = ind_cout) {
        os << "The centralizer is generated by the following element"
           << (number_of_generators() > 1 ? "s" : "") << ":" << EndLine(2);

        os << "───────────" << (number_of_generators() > 1 ? "─" : "")
           << EndLine() << " Generator"
           << (number_of_generators() > 1 ? "s " : " ") << EndLine()
           << "───────────" << (number_of_generators() > 1 ? "─" : "");

        os.Indent(4);

        os << EndLine(1);

        i16 indent =
            (int(std::to_string(number_of_generators() - 1).length()) + 1) / 4 +
            1;
        i16 count = 0;
        for (ConstIterator it = begin(); it != end(); it++) {
            os << count << ":";
            for (size_t _ = 0;
                 _ < 4 * indent - 1 -
                         std::to_string(number_of_generators() - 1).length();
                 _++) {
                os << " ";
            }
            count++;
            os.Indent(4 * indent);
            (*it).print(os);
            os.Indent(-4 * indent);
            os << EndLine();
        }
        os.Indent(-4);
        os << EndLine(1);
    }
};

/**
 * @brief Computes the centralizer of the first element of `uss`.
 *
 * Given an ultra summit set `uss`, with additional structural information
 * stocked in `mins` and `prev`, computes a set of generators of the centralizer
 * of the first braid of `uss` (in the sense that it is the first braid of
 * the first orbit).
 *
 * @tparam F A class representing factors.
 * @param uss The ultra summit set whose first braid's centralizer is to be
 * computed.
 * @param mins A vector that holds, for each `i`, an element that conjugates the
 * base of orbit `prev[i]` to the base of orbit `i`.
 * @param prev A vector that holds integers, such that, for each `i`, `mins[i]`
 * conjugates the base of orbit `prev[i]` to the base of orbit `i`.
 * @return The centralizer of the first braid in `uss`.
 */
template <class F>
Centralizer<BraidTemplate<F>>
centralizer(const ultra_summit::UltraSummitSet<BraidTemplate<F>> &uss,
            const std::vector<F> &mins, const std::vector<i16> &prev) {
    BraidTemplate<F> b = uss.at(0, 0);

    Centralizer<BraidTemplate<F>> centralizer;

    for (size_t orbit_index = 0; orbit_index < uss.number_of_orbits();
         orbit_index++) {
        BraidTemplate<F> d = ultra_summit::tree_path(
                             uss.at(orbit_index, (size_t)0), uss, mins, prev),
                         c = d, b2(b.get_parameter());
        for (size_t shift = 0; shift < uss.orbit_size(orbit_index); shift++) {
            c.right_multiply(
                uss.at(orbit_index, shift).first().delta_conjugate(b.inf()));
        }
        c.right_multiply(!d);

        if (!(c.is_identity())) {
            centralizer.insert(c);
        }

        BraidTemplate<F> orbit_base_rcf = uss.at(orbit_index, (size_t)0);
        orbit_base_rcf.lcf_to_rcf();

        std::vector<F> min = ultra_summit::min_ultra_summit(
            uss.at(orbit_index, (size_t)0), orbit_base_rcf);

        for (typename std::vector<F>::const_iterator it = min.begin();
             it != min.end(); it++) {
            b2 = uss.at(orbit_index, (size_t)0);
            b2.conjugate(*it);
            c = d * (*it) * !ultra_summit::tree_path(b2, uss, mins, prev);

            if (!(c.is_identity())) {
                centralizer.insert(c);
            }
        }
    }

    return centralizer;
}

/**
 * @brief Computes `b`'s centralizer.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose centralizer is to be computed.
 * @return `b`'s centralizer.
 */
template <class F>
Centralizer<BraidTemplate<F>> centralizer(const BraidTemplate<F> &b) {
    std::vector<F> mins;
    std::vector<i16> prev;

    Centralizer<BraidTemplate<F>> centralizer_uss = centralizer(
                                      ultra_summit::ultra_summit_set(b, mins,
                                                                     prev),
                                      mins, prev),
                                  centralizer;

    BraidTemplate<F> c = BraidTemplate<F>(b.get_parameter()), d = c;
    ultra_summit::send_to_ultra_summit(b, c);
    c = !c;

    for (typename Centralizer<BraidTemplate<F>>::ConstIterator it =
             centralizer_uss.begin();
         it != centralizer_uss.end(); it++) {
        d = *it;
        d.conjugate(c);
        centralizer.insert(d);
    }

    return centralizer;
}

} // namespace garcide::centralizer

#endif
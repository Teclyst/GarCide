/**
 * @file super_summit.hpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header (and implementation) file for super summit sets.
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

#ifndef SUPER_SUMMIT
#define SUPER_SUMMIT

#include "garcide/garcide.hpp"
#include <cstddef>
#include <iterator>

/**
 * @brief Namespace for super summit sets.
 */
namespace garcide::super_summit {

/**
 * @brief Computes a super summit conjugate of `b`.
 *
 * This is done through iterated cycling until a repetition is found.
 *
 * @tparam F A class representing factors.
 * @param b The braid of whom a super summit conjugate is computed.
 * @return A super summit conjugate of `b`.
 */
template <class F>
BraidTemplate<F> send_to_super_summit(const BraidTemplate<F> &b) {
    typename F::Parameter n = b.get_parameter();

    i16 k = F(n).lattice_height();

    BraidTemplate<F> b2 = b, b3 = b;

    i16 p = b.inf();
    i16 j = 0;

    while (j <= k) {
        b2.cycling();

        if (b2.inf() == p) {
            j++;
        } else {
            b3 = b2;
            p++;
            j = 0;
        }
    }

    j = 0;
    b2 = b3;
    i16 l = b2.sup();
    while (j <= k) {
        b2.decycling();

        if (b2.sup() == l) {
            j++;
        } else {
            b3 = b2;
            l--;
            j = 0;
        }
    }

    return b3;
}

/**
 * Computes a super summit conjugate of `b`, with a conjugator.
 *
 * @tparam F A class representing factors.
 * @param b The braid of whom a super summit conjugate is computed.
 * @param c A braid that is set to a conjugator sending `b` to what is returned.
 * @return An ultra summit conjugate of `b`.
 */
template <class F>
BraidTemplate<F> send_to_super_summit(const BraidTemplate<F> &b,
                                      BraidTemplate<F> &c) {

    typename F::Parameter n = b.get_parameter();

    i16 k = F(n).lattice_height();

    BraidTemplate<F> b2 = b, b3 = b, c2 = BraidTemplate<F>(n);

    c.identity();

    i16 p = b.inf();
    i16 j = 0;

    while (j <= k) {
        if (b2.canonical_length() == 0) {
            return b2;
        }

        c2.right_multiply(b2.first().delta_conjugate(b2.inf()));
        b2.cycling();

        if (b2.inf() == p) {
            j++;
        } else {
            b3 = b2;
            p++;
            j = 0;
            c.right_multiply(c2);
            c2.identity();
        }
    }

    j = 0;
    b2 = b3;
    i16 l = b2.sup();
    c2.identity();

    while (j <= k) {
        c2.left_multiply(b2.final());
        b2.decycling();

        if (b2.sup() == l) {
            j++;
        } else {
            b3 = b2;
            l--;
            j = 0;
            c.right_divide(c2);
            c2.identity();
        }
    }

    return b3;
}

/**
 * @brief Computes the smallest factor above `f` that conjugates `b` to an
 * element of its summit set.
 *
 * `b` is assumed to be in its summit set.
 *
 * @tparam F A class representing factors.
 * @param b A braid, assumed to be in its summit set.
 * @param f A factor.
 * @return The smallest factor above `f` that conjugates `b` to its summit set.
 */
template <class F> F min_summit(const BraidTemplate<F> &b, const F &f) {
    F r2 = f, r = F(f.get_parameter());
    r.identity();

    BraidTemplate<F> w = b;
    w.set_delta(0);

    while (!r2.is_identity()) {
        r.right_multiply(r2);
        r2 = (w * r).remainder(r.delta_conjugate(b.inf()));
    }

    return r;
}

/**
 * @brief Computes the smallest factor above `f` that conjugates `b` to an
 * element of its super summit set.
 *
 * `b` is assumed to be in its super summit set.
 *
 * @tparam F A class representing factors.
 * @param b A braid, assumed to be in its super summit set.
 * @param b_rcf `b` in RCF.
 * @param f A factor.
 * @return The smallest factor above `f` that conjugates `b` to its super summit
 * set.
 */
template <class F>
F min_super_summit(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
                   const F &f) {
    F r = min_summit(b, f);
    BraidTemplate<F> b2 = b_rcf;
    b2.conjugate_rcf(r);

    while (b2.canonical_length() > b.canonical_length()) {
        r.right_multiply(b2.first());
        b2 = b_rcf;
        b2.conjugate_rcf(r);
    }
    return r;
}

/**
 * @brief Computes the super summit indecomposable conjugators at `b`.
 *
 * `b` is assumed to be in its super summit set.
 *
 * The indecomposable conjugators at `b` are the minimal non-trivial simple
 * factors that conjugate `b` to its super summit set.
 *
 * @tparam F A class representing canonical factors.
 * @param b A braid, assumed to be in its ultra summit set.
 * @param b_rcf `b` in RCF.
 * @return The super summit indecomposable conjugators at `b`.
 */
template <class F>
std::vector<F> min_super_summit(const BraidTemplate<F> &b,
                                const BraidTemplate<F> &b_rcf) {
    F f = F(b.get_parameter());
    std::vector<F> atoms = f.atoms();
    std::vector<F> factors = atoms;

#ifndef USE_PAR

    std::transform(
        atoms.begin(), atoms.end(), factors.begin(),
        [&b, &b_rcf](F &atom) { return min_super_summit(b, b_rcf, atom); });

#else

    std::transform(
        std::execution::par, atoms.begin(), atoms.end(), factors.begin(),
        [&b, &b_rcf](F &atom) { return min_super_summit(b, b_rcf, atom); });

#endif

    std::vector<F> min;

    std::vector<bool> table(atoms.size(), false);
    bool should_be_added;

    for (i16 i = 0; i < int(atoms.size()); i++) {
        f = factors[i];
        should_be_added = true;

        // We check, before adding f, that a divisor of it wasn't added already
        // with some other atom dividing it.
        for (i16 j = 0; j < i && should_be_added; j++) {
            should_be_added = !(table[j] && ((atoms[j] ^ f) == atoms[j]));
        }
        // If that is not the case, we also check if the atom we see is the last
        // that might generate f. This is to avoid duplicates; furthermore, if
        // some strict left divisor of f can be generated by an atom, doing
        // things in this order ensures that it would have been detected by the
        // first loop by the time we try to add f.
        for (i16 j = i + 1; j < int(atoms.size()) && should_be_added; j++) {
            should_be_added = !((atoms[j] ^ f) == atoms[j]);
        }
        if (should_be_added) {
            min.push_back(f);
            table[i] = true;
        }
    }

    return min;
}

/**
 * @brief A class for ultra summit sets.
 *
 * This is wrapper for a `std::unordered_set`.
 *
 * @tparam B A class representing braids.
 */
template <class B> class SuperSummitSet {
  private:
    /**
     * @brief Set of the elements.
     */
    std::unordered_set<B> set;

  public:
    /**
     * @brief Constant iterator type.
     */
    using ConstIterator = typename std::unordered_set<B>::const_iterator;

    /**
     * @brief Constant iterator to the first element of the super summit set.
     *
     * @return An iterator to the first element of `this`.
     */
    inline ConstIterator begin() const { return ConstIterator(set.begin()); }

    /**
     * @brief Constant iterator to the after-last element of the super summit
     * set.
     *
     * @return An iterator to the after-last element of `this`.
     */
    inline ConstIterator end() const { return ConstIterator(set.end()); }

    /**
     * @brief Pushes a braid into the ultra summit set.
     *
     * @param b The braid to be pushed.
     */
    inline void insert(B b) { set.insert(b); }

    /**
     * @brief Membership test.
     *
     * @param b The braid whose membership is tested
     * @return If `b` is in `*this`.
     */
    inline bool mem(const B &b) const { return set.find(b) != set.end(); }

    /**
     * @brief Cardinal of the super summit set.
     *
     * @return The cardinal of `*this`.
     */
    inline i16 card() const { return set.size(); }

    /**
     * @brief Prints the super summit set in output stream `os`.
     *
     * @param os The `IndentedOStream` `*this` is printed in.
     */
    void print(IndentedOStream &os = ind_cout) const {
        os << "There " << (card() > 1 ? "are " : "is ") << card() << " element"
           << (card() > 1 ? "s " : " ") << "in the super summit set."
           << EndLine(2);

        os << "─────" << EndLine() << " Set " << EndLine() << "─────";

        os.Indent(4);

        os << EndLine(1);

        i16 indent = (int(std::to_string(card() - 1).length()) + 1) / 4 + 1;
        i16 count = 0;
        for (ConstIterator it = begin(); it != end(); it++) {
            os << count << ":";
            for (i16 _ = 0;
                 _ < 4 * indent - 1 - int(std::to_string(card() - 1).length());
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

    /**
     * @brief Prints the internal data of the ultra summit set in output stream
     * `os`.
     *
     * Private member `set` is printed.
     *
     * @param os The `IndentedOStream` `*this` is printed in.
     */
    void debug(IndentedOStream &os = ind_cout) const {
        bool is_first = true;
        os << "{   ";
        os.Indent(4);
        os << "set:";
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
};

/**
 * @brief Computes the super summit set of `b`.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose ultra summit set is computed.
 * @return The super summit set of `b`.
 */
template <class F>
SuperSummitSet<BraidTemplate<F>> super_summit_set(const BraidTemplate<F> &b) {
    std::list<BraidTemplate<F>> queue, queue_rcf;
    SuperSummitSet<BraidTemplate<F>> sss;

    BraidTemplate<F> b2 = send_to_super_summit(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.lcf_to_rcf();

    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    sss.insert(b2);

    while (!queue.empty()) {
        std::vector<F> min = min_super_summit(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.conjugate(*itf);

            if (!sss.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.conjugate_rcf(*itf);

                sss.insert(b2);
                queue.push_back(b2);
                queue_rcf.push_back(b2_rcf);
            }
        }

        queue.pop_front();
        queue_rcf.pop_front();
    }

    return sss;
}

/**
 * @brief Checks if two braids are conjugates.
 *
 * This function uses super summit sets.
 *  
 * @tparam F A class representing factors.
 * @param u A braid.
 * @param v Another braid.
 * @return if `u` and `v` are conjugates.
 */
template <class F>
inline bool are_conjugate(const BraidTemplate<F> &u,
                          const BraidTemplate<F> &v) {
    std::unordered_set<BraidTemplate<F>> u_sss = super_summit_set(u);
    return u_sss.mem(send_to_super_summit(v));
}

} // namespace garcide::super_summit

#endif
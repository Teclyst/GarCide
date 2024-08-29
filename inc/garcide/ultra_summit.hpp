/**
 * @file ultra_summit.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header (and implementation) file for ultra summit sets.
 * @version 0.1
 * @date 2024-07-28
 *
 * @copyright Copyright (C) 2024. Distributed under the GNU General Public
 * License, version 3.
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

#ifndef ULTRA_SUMMIT
#define ULTRA_SUMMIT

#include "garcide/super_summit.h"

/**
 * @brief Namespace for ultra summit sets.
 */
namespace garcide::ultra_summit {

/**
 * @brief Exception thrown when a braid that should be ultra summit is not.
 *
 * Exception thrown when a braid that should be ultra summit is not. Used for
 * some functions that may not return when called on non-ultra-summit braids.
 *
 * @tparam B A class representinf braids.
 */
template <class B> struct NotUltraSummit {
    /**
     * @brief The braid that should have been in its ultra summit set.
     */
    B not_ultra_summit;

    /**
     * @brief Construct a new `NotUltraSummit` exception.
     *
     * @param b The braid that should have been in its ultra summit set.
     */
    NotUltraSummit(const B &b) : not_ultra_summit(b) {}
};

/**
 * @brief Computes the trajectory of `b` for cycling.
 *
 * That is, cycles b until a repetition
 * occurs, and then returns the list of the conjugates up to that point.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose trajectory is to be computed.
 * @return The trajectory of `b` for cycling.
 */
template <class F>
std::vector<BraidTemplate<F>> trajectory(BraidTemplate<F> b) {
    std::vector<BraidTemplate<F>> t;
    std::unordered_set<BraidTemplate<F>> t_set;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        b.cycling();
    }

    return t;
}

/**
 * @brief Computes the trajectory of `b` for cycling, in LCF and RCF.
 *
 * That is, cycles b until a repetition
 * occurs, and then returns the list of the conjugates up to that point.
 *
 * This is done in both LCF and RCF simulteanously (as going back and
 * forth is costly).
 *
 * It is assumed that initially `b_rcf` is `b` in RCF.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose trajectory is to be computed.
 * @param b_rcf The braid whose trajectory is to be computed.
 * @param t A vector that the function sets to the trajectory of `b`.
 * @param t_rcf A vector that the function sets to the RCF trajectory of `b`.
 */
template <class F>
void trajectory(BraidTemplate<F> b, BraidTemplate<F> b_rcf,
                std::vector<BraidTemplate<F>> &t,
                std::vector<BraidTemplate<F>> &t_rcf) {
    std::unordered_set<BraidTemplate<F>> t_set;
    t.clear();
    t_rcf.clear();

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_rcf.push_back(b_rcf);
        t_set.insert(b);
        // Cycle in RCF.
        b_rcf.conjugate_rcf(b.initial());
        b.cycling();
    }
}

/**
 * @brief Computes an ultra summit conjugate of `b`.
 *
 * This is done through iterated cycling until a repetition is found.
 *
 * @tparam F A class representing factors.
 * @param b The braid of whom an ultra summit conjugate is computed.
 * @return An ultra summit conjugate of `b`.
 */
template <class F>
BraidTemplate<F> send_to_ultra_summit(const BraidTemplate<F> &b) {
    BraidTemplate<F> b_uss =
        trajectory(super_summit::send_to_super_summit(b)).back();
    b_uss.cycling();
    return b_uss;
}

/**
 * Computes an ultra summit conjugate of `b`, with a conjugator.
 *
 * This is done through iterated cycling until a repetition is found.
 * In the mean time, a conjugator is computed, and `c` is set
 * to it.
 *
 * @tparam F A class representing factors.
 * @param b The braid of whom an ultra summit conjugate is computed.
 * @param c A braid that is set by the function to a conjugator sending `b` to
 * what is returned.
 * @return An ultra summit conjugate of `b`.
 */
template <class F>
BraidTemplate<F> send_to_ultra_summit(const BraidTemplate<F> &b,
                                      BraidTemplate<F> &c) {
    BraidTemplate<F> b_sss = super_summit::send_to_super_summit(b, c);
    std::vector<BraidTemplate<F>> t = trajectory(b_sss);

    BraidTemplate<F> b_uss = BraidTemplate(t.back());
    b_uss.cycling();

    for (typename std::vector<BraidTemplate<F>>::iterator it = t.begin();
         *it != b_uss; it++) {
        c.right_multiply((*it).first().delta_conjugate(b_sss.inf()));
    }

    return b_uss;
}

/**
 * @brief Computes the transport of `f` at `b` for cycling.
 *
 * See Gebhardt, _A New Approach to the Conjugacy Problem in Garside Groups_,
 * 2003, [arXiv:math/0306199 [math.GT]](https://arxiv.org/abs/math/0306199).
 *
 * @tparam F A class representing factors.
 * @param b The braid where a transport is computed.
 * @param f The factor whose transport is computed.
 * @return The transport of `f` at `b`.
 */
template <class F>
BraidTemplate<F> transport(const BraidTemplate<F> &b, const F &f) {
    BraidTemplate<F> b2 = b;
    b2.conjugate(f);
    BraidTemplate<F> b3 = (!BraidTemplate(b.first()) * f) * b2.first();
    return b3.first();
}

/**
 * @brief Computes the iterated transport of `f` at `b` sending `b` to an
 * element in the trajectory of `b.conjugate(f)`.
 *
 * It is assumed that `b` is in its ultra summit set, and that `f` conjugates it
 * to its super summit set.
 *
 * @tparam F A class representing factors.
 * @param b The braid where iterated transports are computed.
 * @param f The factor whose iterated transports are computed.
 * @return The list of the iterated transport of `f` at `b` sending `b` to an
 * element in the trajectory of `b.conjugate(f)`.
 */
template <class F>
std::list<F> transports_sending_to_trajectory(const BraidTemplate<F> &b,
                                              const F &f) {
    std::list<F> ret;
    std::unordered_set<F> ret_set;
    BraidTemplate<F> b1 = BraidTemplate(b),
                     c2 = BraidTemplate<F>(b.get_parameter());
    i16 i, n = 1;
    F f1 = f;

    BraidTemplate<F> c1 = BraidTemplate(b1.first().delta_conjugate(b1.inf()));
    b1.cycling();

    while (b1 != b) {
        c1.right_multiply(BraidTemplate(b1.first().delta_conjugate(b1.inf())));
        b1.cycling();
        n++;
    }

    while (ret_set.find(f1) == ret_set.end()) {
        ret.push_back(f1);
        ret_set.insert(f1);
        b1 = b;
        b1.conjugate(f1);
        c2.identity();
        for (i = 0; i < n; i++) {
            c2.right_multiply(b1.first().delta_conjugate(b1.inf()));
            b1.cycling();
        }

        BraidTemplate<F> b2 = (!c1) * f1 * c2;

        if (b2.inf() == 1) {
            f1.delta();
        } else if (b2.is_identity()) {
            f1.identity();
        } else {
            f1 = b2.first();
        }
    }

    while (ret.front() != f1) {
        ret.pop_front();
    }

    return ret;
}

/**
 * @brief Computes the pullback for cycling of `f` at `b`.
 *
 * See Gebhardt, _A New Approach to the Conjugacy Problem in Garside Groups_,
 * 2003, [arXiv:math/0306199 [math.GT]](https://arxiv.org/abs/math/0306199).
 *
 * @tparam F A class representing factors.
 * @param b The braid where a pullback is computed.
 * @param b_rcf `b`, in RCF.
 * @param f The factor whose pullback is computed.
 * @return The pullback of `f` at `b`.
 */
template <class F>
F pullback(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
           const F &f) {
    F f1 = b.first().delta_conjugate(b.inf() + 1);
    F f2 = f.delta_conjugate();

    BraidTemplate<F> b2 = BraidTemplate(f1) * f2;

    F delta = F(b.get_parameter());
    delta.delta();
    b2.right_multiply(b2.remainder(delta));

    b2.set_delta(b2.inf() - 1);

    F f0 = f;

    if (b2.inf() == 1) {
        f0.delta();
    } else if (b2.is_identity()) {
        f0.identity();
    } else {
        f0 = b2.first();
    }

    F fi = f.delta_conjugate(b.inf());

    for (typename BraidTemplate<F>::ConstFactorItr it = b.cbegin();
         it != b.cend(); it++) {
        if (it != b.cbegin()) {
            fi = fi.left_join(*it) / *it;
        }
    }
    return super_summit::min_super_summit(b, b_rcf, f0.left_join(fi));
}

/**
 * @brief Computes the main pullback for cycling of `f` at `b`.
 *
 * See Gebhardt, _A New Approach to the Conjugacy Problem in Garside Groups_,
 * 2003, [arXiv:math/0306199 [math.GT]](https://arxiv.org/abs/math/0306199).
 *
 * @tparam F A class representing factors.
 * @param b The braid where the main pullback is computed.
 * @param b_rcf `b`, in RCF.
 * @param f The factor whose main pullback is computed.
 * @return The main pullback of `f` at `b`.
 */
template <class F>
F main_pullback(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
                const F &f) {
    std::vector<F> ret;
    std::unordered_map<F, i16> ret_set;

    BraidTemplate<F> b2 = BraidTemplate(b);

    std::vector<BraidTemplate<F>> t, t_rcf;

    trajectory(b, b_rcf, t, t_rcf);

    F f2 = f;
    i16 index = 0;

    while (ret_set.find(f2) == ret_set.end()) {
        ret.push_back(f2);
        ret_set.insert(std::pair(f2, index));
        for (i16 i = int(t.size()) - 1; i >= 0; i--) {
            f2 = pullback(t[i], t_rcf[i], f2);
        }
        index++;
    }
    index = ret_set.at(f2);

    i16 l = ret.size() - index;
    if (index % l == 0) {
        return f2;
    } else {
        return ret[(index / l + 1) * l];
    }
}

/**
 * @brief Computes the smallest factor above `f` that conjugates `b` to an
 * element of its ultra summit set.
 *
 * `b` is assumed to be in its ultra summit set.
 *
 * @tparam F A class representing factors.
 * @param b A braid, assumed to be in its ultra summit set.
 * @param b_rcf `b` in RCF.
 * @param f A factor.
 * @return The smallest factor above `f` that conjugates `b` to its ultra summit
 * set.
 * @exception NotUltraSummit Thrown if `b` is not in its ultra summit set.
 */
template <class F>
F min_ultra_summit(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
                   const F &f) {
    F f2 = super_summit::min_super_summit(b, b_rcf, f);

    std::list<F> ret = transports_sending_to_trajectory(b, f2);

    typename BraidTemplate<F>::FactorItr it;

    for (it = ret.begin(); it != ret.end(); it++) {
        if ((f ^ *it) == f) {
            return *it;
        }
    }

    f2 = main_pullback(b, b_rcf, f);

    ret = transports_sending_to_trajectory(b, f2);

    for (it = ret.begin(); it != ret.end(); it++) {
        if ((f ^ *it) == f) {
            return *it;
        }
    }

    throw NotUltraSummit<BraidTemplate<F>>(b);
}

/**
 * @brief Computes the ultra summit indecomposable conjugators at `b`.
 *
 * `b` is assumed to be in its ultra summit set.
 *
 * The indecomposable conjugators at `b` are the minimal non-trivial simple
 * factors that conjugate `b` to its ultra summit set.
 *
 * @tparam F A class representing factors.
 * @param b A braid, assumed to be in its ultra summit set.
 * @param b_rcf `b` in RCF.
 * @return The ultra summit indecomposable conjugators at `b`.
 */
template <class F>
std::vector<F> min_ultra_summit(const BraidTemplate<F> &b,
                                const BraidTemplate<F> &b_rcf) {
    F f = F(b.get_parameter());
    std::vector<F> atoms = f.atoms();
    std::vector<F> factors = atoms;

#ifndef USE_PAR

    std::transform(
        atoms.begin(), atoms.end(), factors.begin(),
        [&b, &b_rcf](F &atom) { return min_ultra_summit(b, b_rcf, atom); });

#else

    std::transform(
        std::execution::par, atoms.begin(), atoms.end(), factors.begin(),
        [&b, &b_rcf](F &atom) { return min_ultra_summit(b, b_rcf, atom); });

#endif

    std::vector<F> min;

    std::vector<bool> table(atoms.size(), false);
    bool should_be_added;

    for (i16 i = 0; i < int(atoms.size()); i++) {
        f = factors[i];
        should_be_added = true;

        // We check, before adding f, that a divisor of it wasn't added already
        // with some other atom dividing it.
        for (i16 j = 0; (j < i) && should_be_added; j++) {
            should_be_added = !(table[j] && ((atoms[j] ^ f) == atoms[j]));
        }
        // If that is not the case, we also check if the atom we see is the last
        // that might generate f. This is to avoid duplicates; furthermore, if
        // some strict left divisor of f can be generated by an atom, doing
        // things in this order ensures that it would have been detected by the
        // first loop by the time we try to add f.
        for (i16 j = i + 1; (j < int(atoms.size())) && should_be_added;
             j++) {
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
 * @brief Constant iterator class for ultra summit sets.
 *
 * @tparam B A class representing braids.
 */
template <class B> struct UltraSummitConstIterator {

  public:
    /**
     * @brief Iterator category.
     */
    using iterator_category = std::forward_iterator_tag;

    /**
     * @brief Difference type.
     */
    using difference_type = std::ptrdiff_t;

    /**
     * @brief Value type.
     */
    using value_type = B;

    /**
     * @brief Pointer type.
     */
    using pointer = typename std::unordered_map<B, int>::const_iterator;

    /**
     * @brief Reference type.
     */
    using reference = const B &;

  private:
    pointer ptr;

  public:
    /**
     * @brief Constructs a new `UltraSummitConstIterator` from a pointer.
     *
     * @param ptr A pointer to an instance of `B`.
     */
    UltraSummitConstIterator(pointer ptr) : ptr(ptr) {}

    /**
     * @brief Dereference operator.
     *
     * @return A reference to the instance of `B` the iterator points to.
     */
    reference operator*() const { return std::get<0>(*ptr); }

    /**
     * @brief Structure dereference operator.
     *
     * @return The pointer that the iterator corresponds to.
     */
    pointer operator->() { return ptr; }

    /**
     * @brief Prefix incrementation operator.
     *
     * @return A reference to `*this`, after having increased it.
     */
    UltraSummitConstIterator &operator++() {
        ptr++;
        return *this;
    }

    /**
     * @brief Postfix incrementation operator.
     *
     * @return A reference to `*this`, before it was incremented.
     */
    UltraSummitConstIterator operator++(int) {
        UltraSummitConstIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    /**
     * @brief Equality check.
     *
     * @param b Second argument.
     * @return If `*this` and `b` point to the same adress.
     */
    bool operator==(const UltraSummitConstIterator &b) const {
        return ptr == b.ptr;
    }

    /**
     * @brief Unequality check.
     *
     * @param b Second argument.
     * @return If `*this` and `b` do no point to the same adress.
     */
    bool operator!=(const UltraSummitConstIterator &b) const {
        return ptr != b.ptr;
    }
};

/**
 * @brief A class for ultra summit sets.
 *
 * @tparam B A class representing braids.
 */
template <class B> class UltraSummitSet {
  private:
    /**
     * @brief Orbits for cycling.
     */
    std::vector<std::vector<B>> orbits;

    /**
     * @brief Set of the elements.
     *
     * It maps a braid to the index of the orbit it belongs to.
     */
    std::unordered_map<B, i16> set;

  public:
    /**
     * @brief Constant iterator type.
     *
     * Does not iterates through `*this` orbit by orbit.
     */
    using ConstIterator = UltraSummitConstIterator<B>;

    /**
     * @brief Constant iterator to the first element of the ultra summit set.
     *
     * @return An iterator to the first element of `this`.
     */
    inline ConstIterator begin() const { return ConstIterator(set.begin()); }

    /**
     * @brief Constant iterator to the after-last element of the ultra summit
     * set.
     *
     * @return An iterator to the after-last element of `this`.
     */
    inline ConstIterator end() const { return ConstIterator(set.end()); }

    /**
     * @brief Pushes an orbit into the ultra summit set.
     *
     * @param t The orbit to be pushed.
     */
    inline void insert(std::vector<B> t) {
        orbits.push_back(t);
        for (typename std::vector<B>::iterator it = t.begin(); it != t.end();
             it++) {
            set.insert(std::pair(*it, int(orbits.size()) - 1));
        }
    }

    /**
     * @brief Membership test.
     *
     * @param b The braid whose membership is tested
     * @return If `b` is in `*this`.
     */
    inline bool mem(const B &b) const { return set.find(b) != set.end(); }

    /**
     * @brief Access a braid in the ultra summit set with its position.
     *
     * @param orbit_index Index of its orbit.
     * @param shift Position within that orbit.
     * @return The braid at that position
     */
    inline B at(size_t orbit_index, size_t shift) const {
        return orbits[orbit_index][shift];
    }

    /**
     * @brief Finds the orbit of an element of the ultra summit set.
     *
     * @param b The element whose orbit is searched.
     * @return The index of its orbit.
     */
    inline i16 find_orbit(const B &b) const {
        return std::get<1>(*set.find(b));
    }

    /**
     * @brief Number of orbits.
     *
     * @return The number of orbits in `*this`.
     */
    inline size_t number_of_orbits() const { return orbits.size(); }

    /**
     * @brief Cardinal of the ultra summit set.
     *
     * @return The cardinal of `*this`.
     */
    inline size_t card() const { return set.size(); }

    /**
     * @brief Size of a given orbit.
     *
     * @param orbit_index The index of the orbit.
     * @return The size of the orbit.
     */
    inline size_t orbit_size(size_t orbit_index) const {
        return orbits[orbit_index].size();
    }

    /**
     * @brief Size of a given orbit.
     *
     * @param orbit_index The index of the orbit.
     * @return The size of the orbit.
     */
    inline size_t orbit_size(i16 orbit_index) const {
        return orbits[orbit_index].size();
    }

    /**
     * @brief Prints the ultra summit set in output stream `os`.
     *
     * The sizes of the orbits are printed, as well as their contents.
     *
     * @param os The `IndentedOStream` `*this` is printed in.
     */
    void print(IndentedOStream &os = ind_cout) const {

        os << "There " << (card() > 1 ? "are " : "is ") << card() << " element"
           << (card() > 1 ? "s " : " ") << "in the ultra summit set."
           << EndLine(1);

        if (number_of_orbits() > 1) {
            os << "They are split among " << number_of_orbits()
               << " orbits, of respective sizes ";

            for (size_t i = 0; i < number_of_orbits(); i++) {
                os << orbit_size(i)
                   << (i == number_of_orbits() - 1     ? "."
                       : (i == number_of_orbits() - 2) ? " and "
                                                       : ", ");
            }
        } else {
            os << "There is only one orbit.";
        }

        os << EndLine(2);

        for (size_t i = 0; i < number_of_orbits(); i++) {
            std::string str_i = std::to_string(i);
            for (size_t _ = 0; _ < str_i.length() + 8; _++) {
                os << "─";
            }
            os << EndLine() << " orbit " << str_i << EndLine();
            for (size_t _ = 0; _ < str_i.length() + 8; _++) {
                os << "─";
            }
            os.Indent(4);
            os << EndLine(1) << "There " << (orbit_size(i) > 1 ? "are " : "is ")
               << orbit_size(i) << " element"
               << (orbit_size(i) > 1 ? "s " : " ") << "in this orbit."
               << EndLine(1);
            if (orbit_size(i) > 1) {
                os << "They are " << at(i, (size_t)0).rigidity() << "-rigid.";
            } else {
                os << "It is " << at(i, (size_t)0).rigidity() << "-rigid.";
            }
            os << EndLine(1);
            i16 indent =
                (int(std::to_string(orbit_size(i) - 1).length()) + 1) / 4 + 1;
            for (size_t j = 0; j < orbit_size(i); j++) {
                os << j << ":";
                for (size_t _ = 0;
                     _ < 4 * indent - 1 -
                             std::to_string(orbit_size(i) - 1).length();
                     _++) {
                    os << " ";
                }
                os.Indent(4 * indent);
                at(i, j).print(os);
                os.Indent(-4 * indent);
                if (j == orbit_size(i) - 1) {
                    os.Indent(-4);
                } else {
                    os << EndLine();
                }
            }
            os << EndLine(2);
        }
    }

    /**
     * @brief Prints the internal data of the ultra summit set in output stream
     * `os`.
     *
     * Private members `orbits` and `set` are printed.
     *
     * @param os The `IndentedOStream` `*this` is printed in.
     */
    void debug(IndentedOStream &os) const {
        os << "{   ";
        os.Indent(4);
        os << "orbits:";
        os.Indent(4);
        os << EndLine();
        os << "[   ";
        os.Indent(4);
        for (size_t i = 0; i < number_of_orbits(); i++) {
            os << "[   ";
            os.Indent(4);
            for (size_t j = 0; j < number_of_orbits(); j++) {
                at(i, j).debug(os);
                if (j == orbit_size(i) - 1) {
                    os.Indent(-4);
                } else {
                    os << ",";
                }
                os << EndLine();
            }
            os << "]";
            if (i == number_of_orbits() - 1) {
                os.Indent(-4);
            } else {
                os << ",";
            }
            os << EndLine();
        }
        os << "]";
        os.Indent(-4);
        os << EndLine();
        os << "set:";
        os.Indent(4);
        os << EndLine();
        os << "{   ";
        os.Indent(4);
        bool is_first = true;
        for (typename std::unordered_map<B, i16>::const_iterator it =
                 set.begin();
             it != set.end(); it++) {
            if (!is_first) {
                os << "," << EndLine();
            } else {
                is_first = false;
            }
            (*it).first.debug(os);
            os << ": " << (*it).second;
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
 * @brief Computes the ultra summit set of `b`.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose ultra summit set is computed.
 * @return The ultra summit set of `b`.
 */
template <class F>
UltraSummitSet<BraidTemplate<F>> ultra_summit_set(const BraidTemplate<F> &b) {
    UltraSummitSet<BraidTemplate<F>> uss;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    BraidTemplate<F> b2 = send_to_ultra_summit(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.lcf_to_rcf();

    uss.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    while (!queue.empty()) {
        std::vector<F> min = min_ultra_summit(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.conjugate(*itf);

            if (!uss.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.conjugate_rcf(*itf);

                uss.insert(trajectory(b2));
                queue.push_back(b2);
                queue_rcf.push_back(b2_rcf);
            }
        }
        queue.pop_front();
        queue_rcf.pop_front();
    }
    return uss;
}

/**
 * @brief Computes the ultra summit set of `b`, with extra internal structure.
 *
 * `mins` and `prev` are modified to be able to to retrieve conjugators:
 * `prev[i]` is the coordinate of the orbit that is the predecessor the `i`-th
 * one in the graph BFS of its ultra summit set, and `mins[i]` the corresponding
 * conjugator.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose ultra summit set is computed.
 * @param mins A vector that is set to contain, for each `i`, an element that
 * conjugates the base of orbit `prev[i]` to the base of orbit `i`.
 * @param prev A vector that is set to contain integers, such that, for each
 * `i`, `mins[i]` conjugates the base of orbit `prev[i]` to the base of orbit
 * `i`.
 * @return The ultra summit set of `b`.
 */
template <class F>
UltraSummitSet<BraidTemplate<F>> ultra_summit_set(const BraidTemplate<F> &b,
                                                  std::vector<F> &mins,
                                                  std::vector<i16> &prev) {
    UltraSummitSet<BraidTemplate<F>> uss;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    i16 current = 0;
    mins.clear();
    prev.clear();
    mins.push_back(F(b.get_parameter()));
    mins[0].identity();
    prev.push_back(0);

    BraidTemplate<F> b2 = send_to_ultra_summit(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.lcf_to_rcf();

    uss.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    while (!queue.empty()) {
        std::vector<F> min = min_ultra_summit(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.conjugate(*itf);

            if (!uss.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.conjugate_rcf(*itf);

                uss.insert(trajectory(b2));
                queue.push_back(b2);
                queue_rcf.push_back(b2_rcf);

                mins.push_back(*itf);
                prev.push_back(current);
            }
        }
        queue.pop_front();
        queue_rcf.pop_front();

        current++;
    }
    return uss;
}

/**
 * @brief Computes a conjugator from the first element of `uss` to `b`.
 *
 * It is assumed that `b` is an element of `uss`.
 *
 * @tparam F A class representing factors.
 * @param b An element of `uss`.
 * @param uss An ultra summit set.
 * @param mins A vector that holds, for each `i`, an element that conjugates the
 * base of orbit `prev[i]` to the base of orbit `i`.
 * @param prev A vector that holds integers, such that, for each `i`, `mins[i]`
 * conjugates the base of orbit `prev[i]` to the base of orbit `i`.
 * @return A conjugator from the first element of `uss` to `b`.
 */
template <class F>
BraidTemplate<F> tree_path(const BraidTemplate<F> &b,
                           const UltraSummitSet<BraidTemplate<F>> &uss,
                           const std::vector<F> &mins,
                           const std::vector<i16> &prev) {
    BraidTemplate<F> c = BraidTemplate<F>(b.get_parameter());

    if (b.canonical_length() == 0) {
        return c;
    }

    size_t current = (size_t)uss.find_orbit(b);

    for (size_t shift = 0; shift < uss.orbit_size(current); shift++) {
        c.right_multiply(uss.at(current, shift).initial());
    }

    while (current != 0) {
        c.left_multiply(mins[current]);
        current = prev[current];
    }

    return c;
}

/**
 * @brief Checks if two braids are conjugates, and computes a conjugator.
 *
 * `c` is not modified if `b1` and `b2` are not conjugates.
 *
 * This function uses ultra summit sets.
 *
 * @tparam F A class representing factors.
 * @param b1 A braid.
 * @param b2 Another braid.
 * @param c A braid, that is set by the function to the conjugator that takes
 * `b1` to `b2`, if it exists.
 * @return If `b1` and `b2` are conjugates.
 */
template <class F>
bool are_conjugate(const BraidTemplate<F> &b1, const BraidTemplate<F> &b2,
                   BraidTemplate<F> &c) {
    i16 n = b1.get_parameter();
    BraidTemplate<F> c1 = BraidTemplate<F>(n), c2 = BraidTemplate<F>(n);

    BraidTemplate<F> bt1 = send_to_ultra_summit(b1, c1),
                     bt2 = send_to_ultra_summit(b2, c2);

    if (bt1.canonical_length() != bt2.canonical_length() ||
        bt1.sup() != bt2.sup()) {
        return false;
    }

    if (bt1.canonical_length() == 0) {
        c = c1 * !c2;
        return true;
    }

    std::vector<F> mins;
    std::vector<i16> prev;

    UltraSummitSet<BraidTemplate<F>> uss = ultra_summit_set(bt1, mins, prev);

    if (!uss.mem(bt2)) {
        return false;
    }

    c = c1 * tree_path(bt2, uss, mins, prev) * !c2;

    return true;
}
} // namespace garcide::ultra_summit

#endif
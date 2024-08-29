/**
 * @file sliding_circuits.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header (and implementation) file for sliding circuit sets.
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

#ifndef SLIDING_CIRCUITS
#define SLIDING_CIRCUITS

#include "garcide/super_summit.h"

/**
 * @brief Namespace for sliding circuits sets calculations.
 */
namespace garcide::sliding_circuits {

/**
 * @brief Computes the trajectory of `b` for sliding.
 *
 * That is, slides b until a repetition
 * occurs, and then returns the list of the conjugates up to that point.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose trajectory is to be computed.
 * @return The trajectory of `b` for sliding.
 */
template <class F>
std::vector<BraidTemplate<F>> trajectory(BraidTemplate<F> b) {
    std::vector<BraidTemplate<F>> t;
    std::unordered_set<BraidTemplate<F>> t_set;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        b.sliding();
    }

    return t;
}

/**
 * @brief Computes the trajectory of `b` for sliding, with a conjugator and a
 * sliding circuit conjugate.
 *
 * That is, slides b until a repetition
 * occurs, and then returns the list of the conjugates up to that point, with
 * the index of the first sliding circuit conjugate in it. Also returns a
 * conjugator to it.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose trajectory is to be computed.
 * @param c A braid that is set by the function to a conjugator sending `b` to
 * what is returned.
 * @param d An integer that is set by the function to the index of the first
 * sliding circuit conjugate in the trajectory.
 * @return The trajectory of `b` for sliding.
 */
template <class F>
std::vector<BraidTemplate<F>> trajectory(BraidTemplate<F> b,
                                         BraidTemplate<F> &c, i16 &d) {
    std::vector<BraidTemplate<F>> t;
    std::unordered_set<BraidTemplate<F>> t_set;

    c.identity();
    d = 0;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        c.right_multiply(b.preferred_prefix());
        b.sliding();
        d++;
    }

    BraidTemplate<F> b2 = b;
    BraidTemplate<F> c2 = BraidTemplate(b2.preferred_prefix());
    b2.sliding();
    d--;
    while (b2 != b) {
        c2.right_multiply(b2.preferred_prefix());
        b2.sliding();
        d--;
    }
    c.right_divide(c2);

    return t;
}

/**
 * @brief Computes a sliding circuit conjugate of `b`.
 *
 * This is done through iterated sliding until a repetition is found.
 *
 * @tparam F A class representing factors.
 * @param b The braid of whom a sliding circuit conjugate is computed.
 * @return A sliding circuit conjugate of `b`.
 */
template <class F>
BraidTemplate<F> send_to_sliding_circuits(const BraidTemplate<F> &b) {
    BraidTemplate<F> b_sc = trajectory(b).back();
    b_sc.sliding();
    return b_sc;
}

/**
 * Computes a sliding circuit conjugate of `b`, with a conjugator.
 *
 * This is done through iterated sliding until a repetition is found.
 * In the mean time, a conjugator is computed, and `c` is set
 * to it.
 *
 * @tparam F A class representing factors.
 * @param b The braid of whom a sliding circuit conjugate is computed.
 * @param c A braid that is set by the function to a conjugator sending `b` to
 * what is returned.
 * @return A sliding circuit conjugate of `b`.
 */
template <class F>
BraidTemplate<F> send_to_sliding_circuits(const BraidTemplate<F> &b,
                                          BraidTemplate<F> &c) {
    i16 d;
    BraidTemplate<F> b_sc = trajectory(b, c, d).back();
    b_sc.sliding();
    return b_sc;
}

/**
 * @brief Computes the transport of `f` at `b` for sliding.
 *
 * See Gebhardt, González-Meneses, _The cyclic sliding operation in Garside
 * groups_, 2003, [arXiv:0808.1430 [math.GR]](https://arxiv.org/abs/0808.1430).
 *
 * @tparam F A class representing factors.
 * @param b The braid where a transport is computed.
 * @param f The factor whose transport is computed.
 * @return The transport of `f` at `b`.
 */
template <class F> F transport(const BraidTemplate<F> &b, const F &f) {
    BraidTemplate<F> b2 = b;
    b2.conjugate(f);
    BraidTemplate<F> b3 =
        !BraidTemplate(b.preferred_prefix()) * f * b2.preferred_prefix();

    F f2 = F(b.get_parameter());

    if (b3.canonical_length() > 0) {
        f2 = b3.first();
    } else if (b3.inf() == 1) {
        f2.delta();
    } else {
        f2.identity();
    }

    return f2;
}

/**
 * @brief Computes the iterated transport of `f` at `b` sending `b` to an
 * element in the trajectory of `b.conjugate(f)`.
 *
 * It is assumed that `b` is in its sliding circuits set, and that `f`
 * conjugates it to its super summit set.
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
    F g = f;

    BraidTemplate<F> b1 = b;
    i16 i, N = 1;

    b1.sliding();

    while (b1 != b) {
        N++;
        b1.sliding();
    }

    while (ret_set.find(g) == ret_set.end()) {
        ret.push_back(g);
        ret_set.insert(g);

        b1 = b;
        for (i = 0; i < N; i++) {
            g = transport(b1, g);
            b1.sliding();
        }
    }

    while (ret.front() != g) {
        ret.pop_front();
    }

    return ret;
}

/**
 * @brief Computes the pullback for sliding of `f` at `b`.
 *
 * See Gebhardt, González-Meneses, _Solving the Conjugacy Problem in Garside
 * Groups by Cyclic Sliding_, 2003, [arXiv:0809.0948
 * [math.GT]](https://arxiv.org/abs/0809.0948).
 *
 * @tparam F A class representing factors.
 * @param b The braid where a pullback is computed.
 * @param f The factor whose pullback is computed.
 * @return The pullback of `f` at `b`.
 */
template <class F> F pullback(const BraidTemplate<F> &b, const F &f) {
    BraidTemplate<F> b2 = BraidTemplate(b.preferred_prefix());
    b2.right_multiply(f);
    BraidTemplate<F> b3 = b;
    b3.sliding();
    b3.conjugate(f);
    F f2 = b3.preferred_suffix();

    BraidTemplate<F> c = b2.right_meet(f2);

    b2.right_divide(c);

    if (b2.is_identity()) {
        f2.identity();
        return f2;
    } else if (b2.canonical_length() == 0) {
        f2.delta();
        return f2;
    } else {
        return b2.first();
    }
}

/**
 * @brief Computes the main pullback for sliding of `f` at `b`.
 *
 * See Gebhardt, González-Meneses, _Solving the Conjugacy Problem in Garside
 * Groups by Cyclic Sliding_, 2003, [arXiv:0809.0948
 * [math.GT]](https://arxiv.org/abs/0809.0948).
 *
 * @tparam F A class representing factors.
 * @param b The braid where the main pullback is computed.
 * @param f The factor whose main pullback is computed.
 * @return The main pullback of `f` at `b`.
 */
template <class F> F main_pullback(const BraidTemplate<F> &b, const F &f) {
    std::vector<F> ret;
    std::unordered_set<F> ret_set;

    BraidTemplate<F> b2 = b;

    std::vector<BraidTemplate<F>> t = trajectory(b);

    if (f.is_delta()) {
        return f;
    }

    F f2 = f;
    while (ret_set.find(f2) == ret_set.end()) {
        ret.push_back(f2);
        ret_set.insert(f2);

        for (typename std::vector<BraidTemplate<F>>::reverse_iterator itb =
                 t.rbegin();
             itb != t.rend(); itb++) {
            f2 = pullback(*itb, f2);
        }
    }
    return f2;
}

/**
 * @brief Computes the smallest factor above `f` that conjugates `b` to an
 * element of its sliding circuits set.
 *
 * `b` is assumed to be in its sliding circuits set.
 *
 * @tparam F A class representing factors.
 * @param b A braid, assumed to be in its sliding circuits set.
 * @param b_rcf `b` in RCF.
 * @param f A factor.
 * @return The smallest factor above `f` that conjugates `b` to its sliding
 * circuits set.
 */
template <class F>
F min_sliding_circuits(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
                       const F &f) {
    F f2 = super_summit::min_super_summit(b, b_rcf, f);

    std::list<F> ret = transports_sending_to_trajectory(b, f2);
    for (typename std::list<F>::iterator it = ret.begin(); it != ret.end();
         it++) {
        if ((f ^ *it) == f) {
            return *it;
        }
    }

    f2 = main_pullback(b, f);

    ret = transports_sending_to_trajectory(b, f2);

    for (typename std::list<F>::iterator it = ret.begin(); it != ret.end();
         it++) {
        if ((f ^ *it) == f) {
            return *it;
        }
    }

    f2.delta();

    return f2;
}

/**
 * @brief Computes the sliding circuit indecomposable conjugators at `b`.
 *
 * `b` is assumed to be in its sliding circuits set.
 *
 * The indecomposable conjugators at `b` are the minimal non-trivial simple
 * factors that conjugate `b` to its sliding circuits set.
 *
 * @tparam F A class representing factors.
 * @param b A braid, assumed to be in its sliding circuits set.
 * @param b_rcf `b` in RCF.
 * @return The sliding circuit indecomposable conjugators at `b`.
 */
template <class F>
std::vector<F> min_sliding_circuits(const BraidTemplate<F> &b,
                                    const BraidTemplate<F> &b_rcf) {
    F f = F(b.get_parameter());
    std::vector<F> atoms = f.atoms();
    std::vector<F> factors = atoms;

#ifndef USE_PAR

    std::transform(
        atoms.begin(), atoms.end(), factors.begin(),
        [&b, &b_rcf](F &atom) { return min_sliding_circuits(b, b_rcf, atom); });

#else

    std::transform(
        std::execution::par, atoms.begin(), atoms.end(), factors.begin(),
        [&b, &b_rcf](F &atom) { return min_sliding_circuits(b, b_rcf, atom); });

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
 * @brief Constant iterator class for sliding circuits sets.
 *
 * @tparam B A class representing braids.
 */
template <class B> struct SlidingCircuitsConstIterator {

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
     * @brief Constructs a new `SlidingCircuitsConstIterator` from a pointer.
     *
     * @param ptr A pointer to an instance of `B`.
     */
    SlidingCircuitsConstIterator(pointer ptr) : ptr(ptr) {}

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
    SlidingCircuitsConstIterator &operator++() {
        ptr++;
        return *this;
    }

    /**
     * @brief Postfix incrementation operator.
     *
     * @return A reference to `*this`, before it was incremented.
     */
    SlidingCircuitsConstIterator operator++(int) {
        SlidingCircuitsConstIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    /**
     * @brief Equality check.
     *
     * @param b Second argument.
     * @return If `*this` and `b` point to the same adress.
     */
    bool operator==(const SlidingCircuitsConstIterator &b) const {
        return ptr == b.ptr;
    }

    /**
     * @brief Unequality check.
     *
     * @param b Second argument.
     * @return If `*this` and `b` do no point to the same adress.
     */
    bool operator!=(const SlidingCircuitsConstIterator &b) const {
        return ptr != b.ptr;
    }
};

/**
 * @brief A class for sliding circuits sets.
 *
 * @tparam B A class representing braids.
 */
template <class B> class SlidingCircuitsSet {
  private:
    /**
     * @brief Circuits for sliding.
     */
    std::vector<std::vector<B>> circuits;

    /**
     * @brief Set of the elements.
     *
     * It maps a braid to the index of the circuit it belongs to.
     */
    std::unordered_map<B, i16> set;

  public:
    /**
     * @brief Constant iterator type.
     *
     * Does not iterates through `*this` circuit by circuit.
     */
    using ConstIterator = SlidingCircuitsConstIterator<B>;

    /**
     * @brief Constant iterator to the first element of the sliding circuits
     * set.
     *
     * @return An iterator to the first element of `this`.
     */
    inline ConstIterator begin() const { return ConstIterator(set.begin()); }

    /**
     * @brief Constant iterator to the after-last element of the sliding
     * circuits set.
     *
     * @return An iterator to the after-last element of `this`.
     */
    inline ConstIterator end() const { return ConstIterator(set.end()); }

    /**
     * @brief Pushes a circuit into the sliding circuits set.
     *
     * @param t The circuit to be pushed.
     */
    inline void insert(std::vector<B> t) {
        circuits.push_back(t);
        for (typename std::vector<B>::iterator it = t.begin(); it != t.end();
             it++) {
            set.insert(std::pair(*it, int(circuits.size()) - 1));
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
     * @brief Access a braid in the sliding circuits set with its position.
     *
     * @param circuit_index Index of its orbit.
     * @param shift Position within that orbit.
     * @return The braid at that position
     */
    inline B at(size_t circuit_index, size_t shift) const {
        return circuits[circuit_index][shift];
    }

    /**
     * @brief Finds the circuit of an element of the ultra summit set.
     *
     * @param b The element whose circuit is searched.
     * @return The index of its circuit.
     */
    inline i16 find_circuit(const B &b) const { return set.at(b); }

    /**
     * @brief Number of circuits.
     *
     * @return The number of circuits in `*this`.
     */
    inline size_t number_of_circuits() const { return circuits.size(); }

    /**
     * @brief Cardinal of the sliding circuits set.
     *
     * @return The cardinal of `*this`.
     */
    inline size_t card() const { return set.size(); }

    /**
     * @brief Size of a given circuit.
     *
     * @param orbit_index The index of the circuit.
     * @return The size of the circuit.
     */
    inline size_t circuit_size(size_t orbit_index) const {
        return circuits[orbit_index].size();
    }

    /**
     * @brief Size of a given circuit.
     *
     * @param orbit_index The index of the circuit.
     * @return The size of the circuit.
     */
    inline size_t circuit_size(i16 orbit_index) const {
        return circuits[orbit_index].size();
    }

    /**
     * @brief Prints the sliding circuits set in output stream `os`.
     *
     * The sizes of the circuits are printed, as well as their contents.
     *
     * @param os The `IndentedOStream` `*this` is printed in.
     */
    void print(IndentedOStream &os = ind_cout) const {

        os << "There " << (card() > 1 ? "are " : "is ") << card() << " element"
           << (card() > 1 ? "s " : " ") << "in the sliding circuit set."
           << EndLine(1);

        if (number_of_circuits() > 1) {
            os << "They are split among " << number_of_circuits()
               << " circuits, of respective sizes ";

            for (size_t i = 0; i < number_of_circuits(); i++) {
                os << circuit_size(i)
                   << (i == circuits.size() - 1     ? "."
                       : (i == circuits.size() - 2) ? " and "
                                                    : ", ");
            }
        } else {
            os << "There is only one cicuit.";
        }

        os << EndLine(2);

        for (size_t i = 0; i < number_of_circuits(); i++) {
            std::string str_i = std::to_string(i);
            for (size_t _ = 0; _ < str_i.length() + 10; _++) {
                os << "─";
            }
            os << EndLine() << " Circuit " << str_i << EndLine();
            for (size_t _ = 0; _ < str_i.length() + 10; _++) {
                os << "─";
            }
            os.Indent(4);
            os << EndLine(1) << "There "
               << (circuit_size(i) > 1 ? "are " : "is ") << circuit_size(i)
               << " element" << (circuit_size(i) > 1 ? "s " : " ")
               << "in this circuit." << EndLine(1);
            i16 indent =
                (int(std::to_string(circuit_size(i) - 1).length()) + 1) / 4 + 1;
            for (size_t j = 0; j < circuit_size(i); j++) {
                os << j << ":";
                for (size_t _ = 0;
                     _ < 4 * indent - 1 -
                             std::to_string(circuit_size(i) - 1).length();
                     _++) {
                    os << " ";
                }
                os.Indent(4 * indent);
                at(i, j).print(os);
                os.Indent(-4 * indent);
                if (j == circuit_size(i) - 1) {
                    os.Indent(-4);
                } else {
                    os << EndLine();
                }
            }
            os << EndLine(2);
        }
    }

    /**
     * @brief Prints the internal data of the sliding circuits set in output
     * stream `os`.
     *
     * Private members `circuits` and `set` are printed.
     *
     * @param os The `IndentedOStream` `*this` is printed in.
     */
    void debug(IndentedOStream &os) const {
        os << "{   ";
        os.Indent(4);
        os << "circuits:";
        os.Indent(4);
        os << EndLine();
        os << "[   ";
        os.Indent(4);
        for (size_t i = 0; i < number_of_circuits(); i++) {
            os << "[   ";
            os.Indent(4);
            for (size_t j = 0; j < circuit_size(i); j++) {
                at(i, j).debug(os);
                if (j == circuit_size(i) - 1) {
                    os.Indent(-4);
                } else {
                    os << ",";
                }
                os << EndLine();
            }
            os << "]";
            if (i == number_of_circuits() - 1) {
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
 * @brief Computes the sliding circuits set of `b`.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose sliding circuits set is computed.
 * @return The sliding circuits set of `b`.
 */
template <class F>
SlidingCircuitsSet<BraidTemplate<F>>
sliding_circuits_set(const BraidTemplate<F> &b) {
    SlidingCircuitsSet<BraidTemplate<F>> scs;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    BraidTemplate<F> b2 = send_to_sliding_circuits(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.lcf_to_rcf();

    scs.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    F delta = F(b.get_parameter());
    delta.delta();

    b2.conjugate(delta);

    if (!scs.mem(b2)) {
        b2_rcf.conjugate_rcf(delta);

        scs.insert(trajectory(b2));
        queue.push_back(b2);
        queue_rcf.push_back(b2_rcf);
    }

    while (!queue.empty()) {
        std::vector<F> min =
            min_sliding_circuits(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();

            b2.conjugate(*itf);

            if (!scs.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.conjugate_rcf(*itf);

                scs.insert(trajectory(b2));
                queue.push_back(b2);
                queue_rcf.push_back(b2_rcf);

                b2.conjugate(delta);
                if (!scs.mem(b2)) {
                    b2_rcf.conjugate_rcf(delta);

                    scs.insert(trajectory(b2));
                    queue.push_back(b2);
                    queue_rcf.push_back(b2_rcf);
                }
            }
        }
        queue.pop_front();
        queue_rcf.pop_front();
    }
    return scs;
}

/**
 * @brief Computes the sliding circuits set of `b`, with extra internal
 * structure.
 *
 * `mins` and `prev` are modified to be able to to retrieve conjugators:
 * `prev[i]` is the coordinate of the circuit that is the predecessor the `i`-th
 * one in the graph BFS of its sliding circuits set, and `mins[i]` the
 * corresponding conjugator.
 *
 * @tparam F A class representing factors.
 * @param b The braid whose sliding circuits set is computed.
 * @param mins A vector that is set to contain, for each `i`, an element that
 * conjugates the base of circuit `prev[i]` to the base of circuit `i`.
 * @param prev A vector that is set to contain integers, such that, for each
 * `i`, `mins[i]` conjugates the base of circuit `prev[i]` to the base of
 * circuit `i`.
 * @return The sliding circuits set of `b`.
 */
template <class F>
SlidingCircuitsSet<BraidTemplate<F>>
sliding_circuits_set(const BraidTemplate<F> &b, std::vector<F> &mins,
                     std::vector<i16> &prev) {
    SlidingCircuitsSet<BraidTemplate<F>> scs;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    i16 current = 0;
    mins.clear();
    prev.clear();
    mins.push_back(F(b.get_parameter()));
    mins[0].identity();
    prev.push_back(0);

    BraidTemplate<F> b2 = send_to_sliding_circuits(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.lcf_to_rcf();

    scs.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    while (!queue.empty()) {
        std::vector<F> min =
            min_sliding_circuits(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.conjugate(*itf);

            if (!scs.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.conjugate_rcf(*itf);

                scs.insert(trajectory(b2));
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
    return scs;
}

/**
 * @brief Computes a conjugator from the first element of `scs` to `b`.
 *
 * It is assumed that `b` is an element of `scs`.
 *
 * @tparam F A class representing factors.
 * @param b An element of `scs`.
 * @param scs A sliding circuits set.
 * @param mins A vector that is set to contain, for each `i`, an element that
 * conjugates the base of circuit `prev[i]` to the base of circuit `i`.
 * @param prev A vector that is set to contain integers, such that, for each
 * `i`, `mins[i]` conjugates the base of circuit `prev[i]` to the base of
 * circuit `i`.
 * @return A conjugator from the first element of `scs` to `b`.
 */
template <class F>
BraidTemplate<F> tree_path(const BraidTemplate<F> &b,
                           const SlidingCircuitsSet<BraidTemplate<F>> &scs,
                           const std::vector<F> &mins,
                           const std::vector<i16> &prev) {
    BraidTemplate<F> c = BraidTemplate<F>(b.get_parameter());

    if (b.canonical_length() == 0) {
        return c;
    }

    i16 current = scs.find_circuit(b);

    for (size_t shift = 0; shift < scs.circuit_size(current); shift++) {
        c.right_multiply(scs.at(current, shift).preferred_prefix());
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
 * This function uses sliding circuits sets.
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
    typename F::Parameter n = b1.get_parameter();
    BraidTemplate<F> c1 = BraidTemplate<F>(n), c2 = BraidTemplate<F>(n);

    BraidTemplate<F> bt1 = send_to_sliding_circuits(b1, c1),
                     bt2 = send_to_sliding_circuits(b2, c2);

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

    SlidingCircuitsSet<BraidTemplate<F>> scs =
        sliding_circuits_set(bt1, mins, prev);

    if (!scs.mem(bt2)) {
        return false;
    }

    c = c1 * tree_path(bt2, scs, mins, prev) * !c2;

    return true;
}

} // namespace garcide::sliding_circuits

#endif
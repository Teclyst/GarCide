#include "cgarside.h"
#include "garsidesss.h"

namespace cgarside::sliding_circuits {

template <class F>
std::vector<BraidTemplate<F>> trajectory(BraidTemplate<F> b) {
    std::vector<BraidTemplate<F>> t;
    std::unordered_set<BraidTemplate<F>> t_set;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        b.Sliding();
    }

    return t;
}

template <class F>
std::vector<BraidTemplate<F>> trajectory(BraidTemplate<F> b,
                                         BraidTemplate<F> &c, sint16 &d) {
    std::vector<BraidTemplate<F>> t;
    std::unordered_set<BraidTemplate<F>> t_set;

    c.Identity();
    d = 0;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        c.RightProduct(b.PreferredPrefix());
        b.Sliding();
        d++;
    }

    BraidTemplate<F> b2 = b;
    BraidTemplate<F> c2 = BraidTemplate(b2.PreferredPrefix());
    b2.Sliding();
    d--;
    while (b2 != b) {
        c2.RightProduct(b2.PreferredPrefix());
        b2.Sliding();
        d--;
    }
    c.RightDivide(c2);

    return t;
}

template <class F> BraidTemplate<F> send_to_sliding_circuits(const BraidTemplate<F> &b) {
    BraidTemplate<F> b_sc = trajectory(b).back();
    b_sc.Sliding();
    return b_sc;
}

template <class F>
BraidTemplate<F> send_to_sliding_circuits(const BraidTemplate<F> &b, BraidTemplate<F> &c) {
    sint16 d;
    BraidTemplate<F> b_sc = trajectory(b, c, d).back();
    b_sc.Sliding();
    return b_sc;
}

template <class F> F transport(const BraidTemplate<F> &b, const F &f) {
    BraidTemplate<F> b2 = b;
    b2.Conjugate(f);
    BraidTemplate<F> b3 =
        !BraidTemplate(b.PreferredPrefix()) * f * b2.PreferredPrefix();

    F f2 = F(b.GetParameter());

    if (b3.CanonicalLength() > 0) {
        f2 = b3.FactorList.front();
    } else if (b3.Delta == 1) {
        f2.Delta();
    } else {
        f2.Identity();
    }

    return f2;
}

template <class F> std::list<F> transports_sending_to_trajectory(const BraidTemplate<F> &b, const F &f) {
    std::list<F> ret;
    std::unordered_set<F> ret_set;
    F g = f;

    BraidTemplate<F> b1 = b;
    sint16 i, N = 1;

    b1.Sliding();

    while (b1 != b) {
        N++;
        b1.Sliding();
    }

    while (ret_set.find(g) == ret_set.end()) {
        ret.push_back(g);
        ret_set.insert(g);

        b1 = b;
        for (i = 0; i < N; i++) {
            g = transport(b1, g);
            b1.Sliding();
        }
    }

    while (ret.front() != g) {
        ret.pop_front();
    }

    return ret;
}

template <class F> F pullback(const BraidTemplate<F> &b, const F &f) {
    BraidTemplate<F> b2 = BraidTemplate(b.PreferredPrefix());
    b2.RightProduct(f);
    BraidTemplate<F> b3 = b;
    b3.Sliding();
    b3.Conjugate(f);
    F f2 = b3.PreferredSuffix();

    BraidTemplate<F> c = b2.RightMeet(f2);

    b2.RightDivide(c);

    if (b2.IsIdentity()) {
        f2.Identity();
        return f2;
    } else if (b2.CanonicalLength() == 0) {
        f2.Delta();
        return f2;
    } else {
        return b2.FactorList.front();
    }
}

template <class F> F main_pullback(const BraidTemplate<F> &b, const F &f) {
    std::vector<F> ret;
    std::unordered_set<F> ret_set;

    BraidTemplate<F> b2 = b;

    std::vector<BraidTemplate<F>> t = trajectory(b);

    if (f.IsDelta()) {
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

template <class F>
F min_sliding_circuits(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf, const F &f) {
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

    f2.Delta();

    return f2;
}

template <class F>
std::vector<F> min_sliding_circuits(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf) {
    F f = F(b.GetParameter());
    std::vector<F> atoms = f.Atoms();
    std::vector<F> factors = atoms;

#ifndef USE_PAR

    std::transform(atoms.begin(), atoms.end(), factors.begin(),
                   [&b, &b_rcf](F &atom) { return min_sliding_circuits(b, b_rcf, atom); });

#else

    std::transform(std::execution::par, atoms.begin(), atoms.end(),
                   factors.begin(),
                   [&b, &b_rcf](F &atom) { return min_sliding_circuits(b, b_rcf, atom); });

#endif

    std::vector<F> min;

    std::vector<bool> table(atoms.size(), false);
    bool should_be_added;

    for (sint16 i = 0; i < int(atoms.size()); i++) {
        f = factors[i];
        should_be_added = true;

        // We check, before adding f, that a divisor of it wasn't added already
        // with some other atom dividing it.
        for (sint16 j = 0; j < i && should_be_added; j++) {
            should_be_added = !(table[j] && ((atoms[j] ^ f) == atoms[j]));
        }
        // If that is not the case, we also check if the atom we see is the last
        // that might generate f. This is to avoid duplicates; furthermore, if
        // some strict left divisor of f can be generated by an atom, doing
        // things in this order ensures that it would have been detected by the
        // first loop by the time we try to add f.
        for (sint16 j = i + 1; (j < int(atoms.size())) && should_be_added;
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

template <class B> struct SCSConstIterator {

  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = B;
    using pointer = typename std::unordered_map<B, int>::const_iterator;
    using reference = const B &;

  private:
    pointer ptr;

  public:
    SCSConstIterator(pointer ptr) : ptr(ptr) {}

    reference operator*() const { return std::get<0>(*ptr); }
    pointer operator->() { return ptr; }

    // Prefix increment
    SCSConstIterator &operator++() {
        ptr++;
        return *this;
    }

    // Postfix increment
    SCSConstIterator operator++(int) {
        SCSConstIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator==(const SCSConstIterator &b) const { return ptr == b.ptr; }
    bool operator!=(const SCSConstIterator &b) const { return ptr != b.ptr; }
};

// A SCS is stored as both an union of (disjoint) circuits, and a set.
// The set is actually a map: for each key it stores the circuit it belongs to
// (as an index referring to circuits). Both are built concurrently; the set
// part is used to speed up membership tests. Up to names, this is exactly the
// same data structure as `ultra_summit_set::UltraSummitSet`.
template <class B> class SlidingCircuitsSet {
  public:
    std::vector<std::vector<B>> circuits;
    std::unordered_map<B, sint16> set;
    using ConstIterator = SCSConstIterator<B>;

    inline ConstIterator begin() const { return ConstIterator(set.begin()); }

    inline ConstIterator end() const { return ConstIterator(set.end()); }

    // Adds a trajectory to the SCS.
    // Linear in the trajectory's length.
    inline void insert(std::vector<B> t) {
        circuits.push_back(t);
        for (typename std::vector<B>::iterator it = t.begin(); it != t.end();
             it++) {
            set.insert(std::pair(*it, int(circuits.size()) - 1));
        }
    }

    // Checks membership.
    inline bool mem(const B &b) const { return set.find(b) != set.end(); }

    // Finds b's circuit.
    inline sint16 circuit(const B &b) const { return set.at(b); }

    inline size_t number_of_circuits() const { return circuits.size(); }

    inline size_t card() const { return set.size(); }

    inline std::vector<size_t> circuit_sizes() const {
        std::vector<size_t> sizes;
        for (size_t i = 0; i < circuits.size(); i++) {
            sizes.push_back(circuits[i].size());
        }
        return sizes;
    }

    void print(IndentedOStream &os = ind_cout) const {

        std::vector<size_t> sizes = circuit_sizes();
        os << "There " << (card() > 1 ? "are " : "is ") << card() << " element"
           << (card() > 1 ? "s " : " ") << "in the sliding circuit set."
           << EndLine(1);

        if (number_of_circuits() > 1) {
            os << "They are split among " << number_of_circuits()
               << " circuits, of respective sizes ";

            for (sint16 i = 0; i < int(circuits.size()); i++) {
                os << sizes[i]
                   << (i == int(circuits.size()) - 1     ? "."
                       : (i == int(circuits.size()) - 2) ? " and "
                                                         : ", ");
            }
        } else {
            os << "There is only one cicuit.";
        }

        os << EndLine(2);

        for (sint16 i = 0; i < int(circuits.size()); i++) {
            std::string str_i = std::to_string(i);
            for (size_t _ = 0; _ < str_i.length() + 10; _++) {
                os << "─";
            }
            os << EndLine() << " Circuit " << str_i << EndLine();
            for (size_t _ = 0; _ < str_i.length() + 10; _++) {
                os << "─";
            }
            os.Indent(4);
            os << EndLine(1) << "There " << (sizes[i] > 1 ? "are " : "is ")
               << sizes[i] << " element" << (sizes[i] > 1 ? "s " : " ")
               << "in this circuit." << EndLine(1);
            sint16 indent =
                (int(std::to_string(circuits[i].size() - 1).length()) + 1) / 4 +
                1;
            for (sint16 j = 0; j < int(circuits[i].size()); j++) {
                os << j << ":";
                for (sint16 _ = 0;
                     _ <
                     4 * indent - 1 -
                         int(std::to_string(circuits[i].size() - 1).length());
                     _++) {
                    os << " ";
                }
                os.Indent(4 * indent);
                circuits[i][j].Print(os);
                os.Indent(-4 * indent);
                if (j == int(circuits[i].size()) - 1) {
                    os.Indent(-4);
                } else {
                    os << EndLine();
                }
            }
            os << EndLine(2);
        }
    }

    void debug(IndentedOStream &os) const {
        os << "{   ";
        os.Indent(4);
        os << "circuits:";
        os.Indent(4);
        os << EndLine();
        os << "[   ";
        os.Indent(4);
        for (sint16 i = 0; i < int(circuits.size()); i++) {
            os << "[   ";
            os.Indent(4);
            for (sint16 j = 0; j < int(circuits[i].size()); j++) {
                circuits[i][j].Debug(os);
                if (j == int(circuits[i].size()) - 1) {
                    os.Indent(-4);
                } else {
                    os << ",";
                }
                os << EndLine();
            }
            os << "]";
            if (i == int(circuits.size()) - 1) {
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
        for (typename std::unordered_map<B, sint16>::const_iterator it =
                 set.begin();
             it != set.end(); it++) {
            if (!is_first) {
                os << "," << EndLine();
            } else {
                is_first = false;
            }
            (*it).first.Debug(os);
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

template <class F>
SlidingCircuitsSet<BraidTemplate<F>> sliding_circuits_set(const BraidTemplate<F> &b) {
    SlidingCircuitsSet<BraidTemplate<F>> scs;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    BraidTemplate<F> b2 = send_to_sliding_circuits(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.MakeRCFFromLCF();

    scs.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    F delta = F(b.GetParameter());
    delta.Delta();

    b2.Conjugate(delta);

    if (!scs.mem(b2)) {
        b2_rcf.ConjugateRCF(delta);

        scs.insert(trajectory(b2));
        queue.push_back(b2);
        queue_rcf.push_back(b2_rcf);
    }

    while (!queue.empty()) {
        std::vector<F> min = min_sliding_circuits(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();

            b2.Conjugate(*itf);

            if (!scs.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.ConjugateRCF(*itf);

                scs.insert(trajectory(b2));
                queue.push_back(b2);
                queue_rcf.push_back(b2_rcf);

                b2.Conjugate(delta);
                if (!scs.mem(b2)) {
                    b2_rcf.ConjugateRCF(delta);

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

template <class F>
SlidingCircuitsSet<BraidTemplate<F>> sliding_circuits_set(const BraidTemplate<F> &b,
                                        std::vector<F> &mins,
                                        std::vector<sint16> &prev) {
    SlidingCircuitsSet<BraidTemplate<F>> scs;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    sint16 current = 0;
    mins.clear();
    prev.clear();
    mins.push_back(F(b.GetParameter()));
    mins[0].Identity();
    prev.push_back(0);

    BraidTemplate<F> b2 = send_to_sliding_circuits(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.MakeRCFFromLCF();

    scs.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    while (!queue.empty()) {
        std::vector<F> min = min_sliding_circuits(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!scs.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.ConjugateRCF(*itf);

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

template <class F>
BraidTemplate<F> tree_path(const BraidTemplate<F> &b,
                          const SlidingCircuitsSet<BraidTemplate<F>> &scs,
                          const std::vector<F> &mins,
                          const std::vector<sint16> &prev) {
    BraidTemplate<F> c = BraidTemplate<F>(b.GetParameter());

    if (b.CanonicalLength() == 0) {
        return c;
    }

    sint16 current = scs.circuit(b);

    for (typename std::vector<BraidTemplate<F>>::const_iterator itb =
             scs.circuits[current].begin();
         *itb != b; itb++) {
        c.RightProduct((*itb).PreferredPrefix());
    }

    while (current != 0) {
        c.LeftProduct(mins[current]);
        current = prev[current];
    }

    return c;
}

template <class F>
bool are_conjugate(const BraidTemplate<F> &b1, const BraidTemplate<F> &b2,
                  BraidTemplate<F> &c) {
    typename F::ParameterType n = b1.GetParameter();
    BraidTemplate<F> c1 = BraidTemplate<F>(n), c2 = BraidTemplate<F>(n);

    BraidTemplate<F> bt1 = send_to_sliding_circuits(b1, c1), bt2 = send_to_sliding_circuits(b2, c2);

    if (bt1.CanonicalLength() != bt2.CanonicalLength() ||
        bt1.Sup() != bt2.Sup()) {
        return false;
    }

    if (bt1.CanonicalLength() == 0) {
        c = c1 * !c2;
        return true;
    }

    std::vector<F> mins;
    std::vector<sint16> prev;

    SlidingCircuitsSet<BraidTemplate<F>> scs = sliding_circuits_set(bt1, mins, prev);

    if (!scs.mem(bt2)) {
        return false;
    }

    c = c1 * tree_path(bt2, scs, mins, prev) * !c2;

    return true;
}
} // namespace cgarside::sliding_circuits
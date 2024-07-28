#ifndef ULTRA_SUMMIT
#define ULTRA_SUMMIT

#include "garsidesss.h"

// An exception.
struct NotInUSS {};

namespace cgarside::ultra_summit {

template <class F>
std::vector<BraidTemplate<F>> trajectory(BraidTemplate<F> b) {
    std::vector<BraidTemplate<F>> t;
    std::unordered_set<BraidTemplate<F>> t_set;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        b.Cycling();
    }

    return t;
}

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
        b_rcf.ConjugateRCF(b.Initial());
        b.Cycling();
    }
}

template <class F>
BraidTemplate<F> send_to_ultra_summit(const BraidTemplate<F> &b) {
    BraidTemplate<F> b_uss =
        trajectory(super_summit::send_to_super_summit(b)).back();
    b_uss.Cycling();
    return b_uss;
}

template <class F>
BraidTemplate<F> send_to_ultra_summit(const BraidTemplate<F> &b,
                                      BraidTemplate<F> &c) {
    BraidTemplate<F> b_sss = super_summit::send_to_super_summit(b, c);
    std::vector<BraidTemplate<F>> t = trajectory(b_sss);

    BraidTemplate<F> b_uss = BraidTemplate(t.back());
    b_uss.Cycling();

    for (typename std::vector<BraidTemplate<F>>::iterator it = t.begin();
         *it != b_uss; it++) {
        c.right_multiply((*it).First().delta_conjugate(b_sss.Delta));
    }

    return b_uss;
}

template <class F>
BraidTemplate<F> transport(const BraidTemplate<F> &b, const F &f) {
    BraidTemplate<F> b2 = b;
    b2.Conjugate(f);
    BraidTemplate<F> b3 = (!BraidTemplate(b.First()) * f) * b2.First();
    return b3.First();
}

template <class F>
std::list<F> transports_sending_to_trajectory(const BraidTemplate<F> &b,
                                              const F &f) {
    std::list<F> ret;
    std::unordered_set<F> ret_set;
    BraidTemplate<F> b1 = BraidTemplate(b),
                     c2 = BraidTemplate<F>(b.get_parameter());
    sint16 i, n = 1;
    F f1 = f;

    BraidTemplate<F> c1 = BraidTemplate(b1.First().delta_conjugate(b1.Delta));
    b1.Cycling();

    while (b1 != b) {
        c1.right_multiply(BraidTemplate(b1.First().delta_conjugate(b1.Delta)));
        b1.Cycling();
        n++;
    }

    while (ret_set.find(f1) == ret_set.end()) {
        ret.push_back(f1);
        ret_set.insert(f1);
        b1 = b;
        b1.Conjugate(f1);
        c2.identity();
        for (i = 0; i < n; i++) {
            c2.right_multiply(b1.First().delta_conjugate(b1.Delta));
            b1.Cycling();
        }

        BraidTemplate<F> b2 = (!c1) * f1 * c2;

        if (b2.Delta == 1) {
            f1.delta();
        } else if (b2.is_identity()) {
            f1.identity();
        } else {
            f1 = b2.First();
        }
    }

    while (ret.front() != f1) {
        ret.pop_front();
    }

    return ret;
}

template <class F>
F pullback(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
           const F &f) {
    F f1 = b.First().delta_conjugate(b.Delta + 1);
    F f2 = f.delta_conjugate();

    BraidTemplate<F> b2 = BraidTemplate(f1) * f2;

    F delta = F(b.get_parameter());
    delta.delta();
    b2.right_multiply(b2.Remainder(delta));

    b2.Delta--;

    F f0 = f;

    if (b2.Delta == 1) {
        f0.delta();
    } else if (b2.is_identity()) {
        f0.identity();
    } else {
        f0 = b2.First();
    }

    F fi = f.delta_conjugate(b.Delta);

    for (typename BraidTemplate<F>::ConstFactorItr it = b.FactorList.begin();
         it != b.FactorList.end(); it++) {
        if (it != b.FactorList.begin()) {
            fi = fi.left_join(*it) / *it;
        }
    }
    return super_summit::min_super_summit(b, b_rcf, f0.left_join(fi));
}

template <class F>
F main_pullback(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
                const F &f) {
    std::vector<F> ret;
    std::unordered_map<F, sint16> ret_set;

    BraidTemplate<F> b2 = BraidTemplate(b);

    std::vector<BraidTemplate<F>> t, t_rcf;

    trajectory(b, b_rcf, t, t_rcf);

    F f2 = f;
    sint16 index = 0;

    while (ret_set.find(f2) == ret_set.end()) {
        ret.push_back(f2);
        ret_set.insert(std::pair(f2, index));
        for (sint16 i = int(t.size()) - 1; i >= 0; i--) {
            f2 = pullback(t[i], t_rcf[i], f2);
        }
        index++;
    }
    index = ret_set.at(f2);

    sint16 l = ret.size() - index;
    if (index % l == 0) {
        return f2;
    } else {
        return ret[(index / l + 1) * l];
    }
}

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

    throw NotInUSS();
}

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

    for (sint16 i = 0; i < int(atoms.size()); i++) {
        f = factors[i];
        should_be_added = true;

        // We check, before adding f, that a divisor of it wasn't added already
        // with some other atom dividing it.
        for (sint16 j = 0; (j < i) && should_be_added; j++) {
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

template <class B> struct USSConstIterator {

  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = B;
    using pointer = typename std::unordered_map<B, int>::const_iterator;
    using reference = const B &;

  private:
    pointer ptr;

  public:
    USSConstIterator(pointer ptr) : ptr(ptr) {}

    reference operator*() const { return std::get<0>(*ptr); }
    pointer operator->() { return ptr; }

    // Prefix increment
    USSConstIterator &operator++() {
        ptr++;
        return *this;
    }

    // Postfix increment
    USSConstIterator operator++(int) {
        USSConstIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    bool operator==(const USSConstIterator &b) const { return ptr == b.ptr; }
    bool operator!=(const USSConstIterator &b) const { return ptr != b.ptr; }
};

// An USS is stored as, on one side, an union of (disjoint) trajectories, and on
// the other side a set. The set is actually a map: for each key it stores the
// orbit it belongs to (as an index referring to orbits). Both are built
// concurrently; the set part is used to speed up membership tests.
template <class B> class UltraSummitSet {
  public:
    std::vector<std::vector<B>> orbits;
    std::unordered_map<B, sint16> set;

    using ConstIterator = USSConstIterator<B>;

    inline ConstIterator begin() const { return ConstIterator(set.begin()); }

    inline ConstIterator end() const { return ConstIterator(set.end()); }

    // Adds a trajectory to the USS.
    // Linear in the trajectory's length.
    inline void insert(std::vector<B> t) {
        orbits.push_back(t);
        for (typename std::vector<B>::iterator it = t.begin(); it != t.end();
             it++) {
            set.insert(std::pair(*it, int(orbits.size()) - 1));
        }
    }

    // Checks membership.
    inline bool mem(const B &b) const { return set.find(b) != set.end(); }

    /**
     * @brief Access a braid in the USS with its position.
     *
     * @param orbit_index Index of its orbit.
     * @param shift Position within that orbit.
     * @return B
     */
    inline B at(size_t orbit_index, size_t shift) const {
        return orbits[orbit_index][shift];
    }

    /**
     * @brief Access a braid in the USS with its position.
     *
     * @param orbit_index Index of its orbit.
     * @param shift Position within that orbit.
     * @return B
     */
    inline B at(sint16 orbit_index, sint16 shift) const {
        return orbits[orbit_index][shift];
    }

    // Finds b's orbit.
    inline sint16 find_orbit(const B &b) const {
        return std::get<1>(*set.find(b));
    }

    inline size_t number_of_orbits() const { return orbits.size(); }

    inline size_t card() const { return set.size(); }

    inline size_t orbit_size(size_t orbit_index) const {
        return orbits[orbit_index].size();
    }

    inline size_t orbit_size(sint16 orbit_index) const {
        return orbits[orbit_index].size();
    }

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
                os << "They are " << at(i, (size_t)0).Rigidity() << "-rigid.";
            } else {
                os << "It is " << at(i, (size_t)0).Rigidity() << "-rigid.";
            }
            os << EndLine(1);
            sint16 indent =
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
        for (typename std::unordered_map<B, sint16>::const_iterator it =
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

template <class F>
UltraSummitSet<BraidTemplate<F>> ultra_summit_set(const BraidTemplate<F> &b) {
    UltraSummitSet<BraidTemplate<F>> uss;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    BraidTemplate<F> b2 = send_to_ultra_summit(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.MakeRCFFromLCF();

    uss.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    while (!queue.empty()) {
        std::vector<F> min = min_ultra_summit(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!uss.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.ConjugateRCF(*itf);

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

template <class F>
UltraSummitSet<BraidTemplate<F>> ultra_summit_set(const BraidTemplate<F> &b,
                                                  std::vector<F> &mins,
                                                  std::vector<sint16> &prev) {
    UltraSummitSet<BraidTemplate<F>> uss;
    std::list<BraidTemplate<F>> queue, queue_rcf;

    sint16 current = 0;
    mins.clear();
    prev.clear();
    mins.push_back(F(b.get_parameter()));
    mins[0].identity();
    prev.push_back(0);

    BraidTemplate<F> b2 = send_to_ultra_summit(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.MakeRCFFromLCF();

    uss.insert(trajectory(b2));
    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    while (!queue.empty()) {
        std::vector<F> min = min_ultra_summit(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!uss.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.ConjugateRCF(*itf);

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

template <class F>
BraidTemplate<F> tree_path(const BraidTemplate<F> &b,
                           const UltraSummitSet<BraidTemplate<F>> &uss,
                           const std::vector<F> &mins,
                           const std::vector<sint16> &prev) {
    BraidTemplate<F> c = BraidTemplate<F>(b.get_parameter());

    if (b.CanonicalLength() == 0) {
        return c;
    }

    sint16 current = uss.find_orbit(b);

    for (typename std::vector<BraidTemplate<F>>::const_iterator itb =
             uss.orbits[current].begin();
         *itb != b; itb++) {
        c.right_multiply((*itb).First().delta_conjugate(b.Delta));
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
    sint16 n = b1.get_parameter();
    BraidTemplate<F> c1 = BraidTemplate<F>(n), c2 = BraidTemplate<F>(n);

    BraidTemplate<F> bt1 = send_to_ultra_summit(b1, c1),
                     bt2 = send_to_ultra_summit(b2, c2);

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

    UltraSummitSet<BraidTemplate<F>> uss = ultra_summit_set(bt1, mins, prev);

    if (!uss.mem(bt2)) {
        return false;
    }

    c = c1 * tree_path(bt2, uss, mins, prev) * !c2;

    return true;
}
} // namespace cgarside::ultra_summit

#endif
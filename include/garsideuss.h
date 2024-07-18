#include "garsidesss.h"

// An exception.
struct NotInUSS {};

namespace USS {
using namespace CGarside;

template <class F> std::vector<Braid<F>> Trajectory(Braid<F> b) {
    std::vector<Braid<F>> t;
    std::unordered_set<Braid<F>> t_set;

    while (t_set.find(b) == t_set.end()) {
        t.push_back(b);
        t_set.insert(b);
        b.Cycling();
    }

    return t;
}

template <class F> Braid<F> SendToUSS(const Braid<F> &b) {
    Braid<F> b_uss = Trajectory(SSS::SendToSSS(b)).back();
    b_uss.Cycling();
    return b_uss;
}

template <class F> Braid<F> SendToUSS(const Braid<F> &b, Braid<F> &c) {
    Braid<F> b_sss = SSS::SendToSSS(b, c);
    std::list<Braid<F>> t = Trajectory(b_sss);

    Braid<F> b_uss = Braid(t.back());
    b_uss.cycling();

    for (typename std::vector<Braid<F>>::iterator it = t.begin(); *it != b_uss;
         it++) {
        c.RightProduct((*it).DeltaConjugate(b_sss.Delta));
    }

    return b_uss;
}

template <class F> Braid<F> Transport(const Braid<F> &b, const F &f) {
    Braid<F> b2 = b;
    b2.Conjugate(f);
    Braid<F> b3 = (!Braid(b.FactorList.front()) * f) * b2.FactorList.front();
    return b3.FactorList.front();
}

template <class F> std::list<F> Return(const Braid<F> &b, const F &f) {
    std::list<F> ret;
    std::unordered_set<F> ret_set;
    Braid<F> b1 = Braid(b), c2 = Braid<F>(b.GetParameter());
    sint16 i, n = 1;
    F f1 = F(f);

    Braid<F> c1 = Braid(b1.FactorList.front().DeltaConjugate(b1.Delta));
    b1.Cycling();

    while (b1 != b) {
        c1.RightProduct(Braid(b1.FactorList.front().DeltaConjugate(b1.Delta)));
        b1.Cycling();
        n++;
    }

    while (ret_set.find(f1) == ret_set.end()) {
        ret.push_back(f1);
        ret_set.insert(f1);
        b1.Conjugate(f1);
        c2.Identity();
        for (i = 0; i < n; i++) {
            c2.RightProduct(b1.FactorList.front().DeltaConjugate(b1.Delta));
            b1.Cycling();
        }

        Braid<F> b2 = (!c1) * f1 * c2;

        if (b2.Delta == 1) {
            f1.Delta();
        } else if (b2.IsIdentity()) {
            f1.Identity();
        } else {
            f1 = b2.FactorList.front();
        }
    }

    while (ret.front() != f1) {
        ret.pop_front();
    }

    return ret;
}

template <class F> F Pullback(const Braid<F> &b, const F &f) {
    F f1 = b.FactorList.front().DeltaConjugate(b.Delta + 1);
    F f2 = f.DeltaConjugate();

    Braid<F> b2 = Braid(f1) * f2;

    F delta = F(b.GetParameter());
    delta.Delta();
    b2.RightProduct(b2.Remainder(delta));

    b2.Delta--;

    F f0 = f;

    if (b2.Delta == 1) {
        f0.Delta();
    } else if (b2.IsIdentity()) {
        f0.Identity();
    } else {
        f0 = b2.FactorList.front();
    }

    F fi = f.DeltaConjugate(b.Delta);

    for (typename Braid<F>::ConstFactorItr it = b.FactorList.begin();
         it != b.FactorList.end(); it++) {
        if (it != b.FactorList.begin()) {
            fi = fi.LeftJoin(*it) / *it;
        }
    }
    return SSS::MinSSS(b, f0.LeftJoin(fi));
}

template <class F> F MainPullback(const Braid<F> &b, const F &f) {
    std::vector<F> ret;
    std::unordered_map<F, sint16> ret_set;

    Braid<F> b2 = Braid(b);

    std::vector<Braid<F>> t = Trajectory(b);

    F f2 = F(f);
    sint16 index = 0;

    while (ret_set.find(f2) == ret_set.end()) {
        ret.push_back(f2);
        ret_set.insert(std::pair(f2, index));
        for (typename std::vector<Braid<F>>::reverse_iterator revit =
                 t.rbegin();
             revit != t.rend(); revit++) {
            f2 = Pullback(*revit, f2);
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
F MinUSS(const Braid<F> &b, const Braid<F> &b_rcf, const F &f) {
    F f2 = SSS::MinSSS(b, b_rcf, f);

    std::list<F> ret = Return(b, f2);

    typename Braid<F>::FactorItr it;

    for (it = ret.begin(); it != ret.end(); it++) {
        if ((f ^ *it) == f) {
            return *it;
        }
    }

    f2 = MainPullback(b, f);

    ret = Return(b, f2);

    for (it = ret.begin(); it != ret.end(); it++) {
        if ((f ^ *it) == f) {
            return *it;
        }
    }

    throw NotInUSS();
}

template <class F> std::vector<F> MinUSS(const Braid<F> &b) {
    F f = F(b.GetParameter());
    std::vector<F> atoms = f.Atoms();
    Braid<F> b_rcf = b;
    b_rcf.MakeRCFFromLCF();

    std::vector<F> min;

    bool table[atoms.size()] = {false};
    bool should_be_added;

    for (sint16 i = 0; i < int(atoms.size()); i++) {
        f = MinUSS(b, b_rcf, atoms[i]);
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

// An USS is stored as, on one side, an union of (disjoint) trajectories, and on
// the other side a set. The set is actually a map: for each key it stores the
// orbit it belongs to (as an index referring to Orbits). Both are built
// concurrently; the set part is used to speed up membership tests.
template <class B> class UltraSummitSet {
    std::vector<std::vector<B>> Orbits;
    std::unordered_map<B, sint16> Set;

  public:
    // Adds a trajectory to the USS.
    // Linear in the trajectory's length.
    inline void Insert(std::vector<B> t) {
        Orbits.push_back(t);
        for (typename std::vector<B>::iterator it = t.begin(); it != t.end();
             it++) {
            Set.insert(std::pair(*it, int(Orbits.size()) - 1));
        }
    }

    // Checks membership.
    inline bool Mem(const B &b) const { return Set.find(b) != Set.end(); }

    // Finds b's orbit.
    inline sint16 Orbit(const B &b) const { return *Set.find(b); }

    void Print(std::ostream &os) const {
        for (sint16 i = 0; i < int(Orbits.size()); i++) {
            os << "Orbit " << i << ":" << std::endl;
            for (sint16 j = 0; j < int(Orbits[i].size()); j++) {
                Orbits[i][j].Print(os);
                os << std::endl;
            }
        }
    }

    void Debug(std::ostream &os) const {
        for (sint16 i = 0; i < int(Orbits.size()); i++) {
            os << "Orbit " << i << ":" << std::endl;
            for (sint16 j = 0; j < int(Orbits[i].size()); j++) {
                Orbits[i][j].Debug(os);
                os << std::endl;
            }
        }
    }
};

template <class F> UltraSummitSet<Braid<F>> USS(const Braid<F> &b) {
    UltraSummitSet<Braid<F>> uss;
    std::list<Braid<F>> queue;
    F f = F(b.GetParameter());

    Braid<F> b2 = SendToUSS(b);

    uss.Insert(Trajectory(b2));
    queue.push_back(b2);

    F delta = F(b.GetParameter());
    delta.Delta();

    b2.Conjugate(delta);

    if (!uss.Mem(b2)) {
        uss.Insert(Trajectory(b2));
        queue.push_back(b2);
    }

    while (!queue.empty()) {
        std::vector<F> min = MinUSS(queue.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!uss.Mem(b2)) {
                uss.Insert(Trajectory(b2));
                queue.push_back(b2);
                b2.Conjugate(delta);
                if (!uss.Mem(b2)) {
                    uss.Insert(Trajectory(b2));
                    queue.push_back(b2);
                }
            }
        }
        queue.pop_front();
    }
    return uss;
}

template <class F>
UltraSummitSet<Braid<F>> USS(const Braid<F> &b, std::vector<F> &mins,
                             std::vector<sint16> &prev) {
    UltraSummitSet<Braid<F>> uss;
    std::list<Braid<F>> queue;
    F f = F(b.GetParameter());

    sint16 current = 0;
    mins.clear();
    prev.clear();

    Braid<F> b2 = SendToUSS(b);

    uss.Insert(Trajectory(b2));
    queue.push_back(b2);

    while (!queue.empty()) {
        std::vector<F> min = MinUSS(queue.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!uss.Mem(b2)) {
                uss.Insert(Trajectory(b2));
                queue.push_back(b2);
                mins.push_back(*itf);
                prev.push_back(current);
            }
        }
        queue.pop_front();
        current++;
    }
    return uss;
}

template <class F>
Braid<F> TreePath(const Braid<F> &b, const UltraSummitSet<Braid<F>> &uss,
                  const std::vector<F> &mins, const std::vector<sint16> &prev) {
    Braid<F> c = Braid(b.GetParameter());

    if (b.CanonicalLength() == 0) {
        return c;
    }

    sint16 current = uss.Orbit(b);

    for (typename std::list<Braid<F>>::iterator itb =
             uss.Orbits[current].begin();
         *itb != b; itb++) {
        c.RightProduct((*itb).FactorList.front().DeltaConjugate(b.Delta));
    }

    while (current != 0) {
        c.LeftMultiply(mins[current]);
        current = prev[current];
    }

    return c;
}

template <class F>
bool AreConjugate(const Braid<F> &b1, const Braid<F> &b2, Braid<F> &c) {
    sint16 n = b1.GetParameter();
    Braid<F> c1 = Braid<F>(n), c2 = Braid<F>(n);

    Braid<F> bt1 = SendToUSS(b1, c1), bt2 = SendToUSS(b2, c2);

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

    UltraSummitSet<Braid<F>> uss = USS(bt1, mins, prev);

    if (!uss.Mem(bt2)) {
        return false;
    }

    c = c1 * TreePath(bt2, uss, mins, prev) * !c2;

    return true;
}
} // namespace USS
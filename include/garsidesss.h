#ifndef SUPER
#define SUPER

#include "cgarside.h"

namespace SSS {
using namespace CGarside;
using CGarside::sint16;

template <class F> Braid<F> SendToSSS(const Braid<F> &b) {
    typename F::ParameterType n = b.GetParameter();

    sint16 k = F(n).LatticeHeight();

    Braid<F> b2 = b, b3 = b;

    sint16 p = b.Delta;
    sint16 j = 0;

    while (j <= k) {
        b2.Cycling();

        if (b2.Delta == p) {
            j++;
        } else {
            b3 = b2;
            p++;
            j = 0;
        }
    }

    j = 0;
    sint16 l = b3.Sup();
    while (j <= k) {
        b3.Decycling();

        if (b3.Sup() == l) {
            j++;
        } else {
            b2 = b3;
            l--;
            j = 0;
        }
    }

    return b2;
}

template <class F> Braid<F> SendToSSS(const Braid<F> &b, const Braid<F> &c) {

    typename F::ParameterType n = b.GetParameter();

    sint16 k = F(n).LatticeHeight();

    Braid<F> b2 = Braid(b), b3 = Braid(b), c2 = Braid(n);

    c.Identity();

    sint16 p = b.Delta;
    sint16 j = 0;

    while (j <= k) {
        if (b2.CanonicalLength() == 0) {
            return b2;
        }

        c2.RightProduct(b2.FactorList.front().DeltaConjugate(b2.delta));
        b2.Cycling();

        if (b2.LeftDelta == p) {
            j++;
        } else {
            b3 = b2;
            p++;
            j = 0;
            c.RightProduct(c2);
            c2.Identity();
        }
    }

    j = 0;
    sint16 l = b3.Sup();
    c2.Identity();

    while (j <= k) {
        c2.LeftProduct(b3.Final());
        b3.Decycling();

        if (b3.Sup() == l) {
            j++;
        } else {
            b2 = b3;
            l--;
            j = 0;
            c.RightDivide(c2);
            c2.Identity();
        }
    }

    return b2;
}

template <class F> F MinSS(const Braid<F> &b, const F &f) {
    F r2 = f, r = F(f.GetParameter());
    r.Identity();

    Braid<F> w = b;
    w.Delta = 0;

    while (!r2.IsIdentity()) {
        r.RightProduct(r2);
        r2 = (w * r).Remainder(r.DeltaConjugate(b.Delta));
    }

    return r;
}

template <class F>
F MinSSS(const Braid<F> &b, const Braid<F> &b_rcf, const F &f) {
    F r = MinSS(b, f);
    Braid<F> b2 = b_rcf;
    b2.ConjugateRCF(r);

    while (b2.CanonicalLength() > b.CanonicalLength()) {
        r.RightProduct(b2.FactorList.front());
        b2 = b_rcf;
        b2.ConjugateRCF(r);
    }
    return r;
}

template <class F>
std::vector<F> MinSSS(const Braid<F> &b, const Braid<F> &b_rcf) {
    F f = F(b.GetParameter());
    std::vector<F> atoms = f.Atoms();
    std::vector<F> factors = atoms;

    std::transform(execution_policy, atoms.begin(), atoms.end(),
                   factors.begin(),
                   [&b, &b_rcf](F &atom) { return MinSSS(b, b_rcf, atom); });

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
        for (sint16 j = i + 1; j < int(atoms.size()) && should_be_added; j++) {
            should_be_added = !((atoms[j] ^ f) == atoms[j]);
        }
        if (should_be_added) {
            min.push_back(f);
            table[i] = true;
        }
    }

    return min;
}

// A SuperSummitSet is basically a wrapper for an unordered set.
// The reason for doing it that way being to have function signatures that are
// closer to those for USSs and SCSs. This might change later - it could be
// interesting to instead use a map and save for each element a conjugator.
template <class B> class SuperSummitSet {
  private:
    std::unordered_set<B> Set;

  public:
    inline void Insert(B b) { Set.insert(b); }

    // Checks membership.
    inline bool Mem(const B &b) const { return Set.find(b) != Set.end(); }

    inline sint16 Card() const { return Set.size(); }

    void Print(IndentedOStream &os) const {
        bool is_first = true;
        for (typename std::unordered_set<B>::const_iterator it =
                 Set.begin();
             it != Set.end(); it++) {
            if (!is_first) {
                os << EndLine();
            } else {
                is_first = false;
            }
            (*it).Print(os);
        }
    }

    void Debug(IndentedOStream &os) const {
        bool is_first = true;
        os << "{   ";
        os.Indent(4);
        os << "Set:";
        os.Indent(4);
        os << EndLine();
        os << "{   ";
        os.Indent(4);
        for (typename std::unordered_set<B>::const_iterator it =
                 Set.begin();
             it != Set.end(); it++) {
            if (!is_first) {
                os << "," << EndLine();
            } else {
                is_first = false;
            }
            (*it).Debug(os);
        }
        os.Indent(-4);
        os << EndLine();
        os << "}";
        os.Indent(-8);
        os << EndLine();
        os << "}";
    }
};

template <class F> SuperSummitSet<Braid<F>> SSS(const Braid<F> &b) {
    std::list<Braid<F>> queue, queue_rcf;
    SuperSummitSet<Braid<F>> sss;
    F f = F(b.GetParameter());

    Braid<F> b2 = SendToSSS(b);
    Braid<F> b2_rcf = b2;
    b2_rcf.MakeRCFFromLCF();

    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    sss.Insert(b2);

    while (!queue.empty()) {
        std::vector<F> min = MinSSS(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!sss.Mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.ConjugateRCF(*itf);
                
                sss.Insert(b2);
                queue.push_back(b2);
                queue_rcf.push_back(b2_rcf);
            }
        }

        queue.pop_front();
        queue_rcf.pop_front();
    }

    return sss;
}

template <class F>
inline bool AreConjugate(const Braid<F> &u, const Braid<F> &v) {
    std::unordered_set<Braid<F>> u_sss = SSS(u);
    return u_sss.Mem(SendToSSS(v));
}

} // namespace SSS

#endif
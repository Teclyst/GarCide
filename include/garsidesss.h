#ifndef SUPER
#define SUPER

#include "cgarside.h"
#include <cstddef>
#include <iterator>

namespace cgarside::super_summit {

template <class F>
BraidTemplate<F> send_to_super_summit(const BraidTemplate<F> &b) {
    typename F::Parameter n = b.get_parameter();

    sint16 k = F(n).lattice_height();

    BraidTemplate<F> b2 = b, b3 = b;

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
    b2 = b3;
    sint16 l = b2.Sup();
    while (j <= k) {
        b2.Decycling();

        if (b2.Sup() == l) {
            j++;
        } else {
            b3 = b2;
            l--;
            j = 0;
        }
    }

    return b3;
}

template <class F>
BraidTemplate<F> send_to_super_summit(const BraidTemplate<F> &b,
                                      BraidTemplate<F> &c) {

    typename F::Parameter n = b.get_parameter();

    sint16 k = F(n).lattice_height();

    BraidTemplate<F> b2 = b, b3 = b, c2 = BraidTemplate<F>(n);

    c.identity();

    sint16 p = b.Delta;
    sint16 j = 0;

    while (j <= k) {
        if (b2.CanonicalLength() == 0) {
            return b2;
        }

        c2.right_multiply(b2.First().delta_conjugate(b2.Inf()));
        b2.Cycling();

        if (b2.Delta == p) {
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
    sint16 l = b2.Sup();
    c2.identity();

    while (j <= k) {
        c2.LeftProduct(b2.Final());
        b2.Decycling();

        if (b2.Sup() == l) {
            j++;
        } else {
            b3 = b2;
            l--;
            j = 0;
            c.RightDivide(c2);
            c2.identity();
        }
    }

    return b3;
}

template <class F> F min_summit(const BraidTemplate<F> &b, const F &f) {
    F r2 = f, r = F(f.get_parameter());
    r.identity();

    BraidTemplate<F> w = b;
    w.Delta = 0;

    while (!r2.is_identity()) {
        r.right_multiply(r2);
        r2 = (w * r).Remainder(r.delta_conjugate(b.Delta));
    }

    return r;
}

template <class F>
F min_super_summit(const BraidTemplate<F> &b, const BraidTemplate<F> &b_rcf,
                   const F &f) {
    F r = min_summit(b, f);
    BraidTemplate<F> b2 = b_rcf;
    b2.ConjugateRCF(r);

    while (b2.CanonicalLength() > b.CanonicalLength()) {
        r.right_multiply(b2.FactorList.front());
        b2 = b_rcf;
        b2.ConjugateRCF(r);
    }
    return r;
}

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
    std::unordered_set<B> set;

  public:
    using ConstIterator = typename std::unordered_set<B>::const_iterator;

    inline ConstIterator begin() const { return ConstIterator(set.begin()); }

    inline ConstIterator end() const { return ConstIterator(set.end()); }

    inline void insert(B b) { set.insert(b); }

    // Checks membership.
    inline bool mem(const B &b) const { return set.find(b) != set.end(); }

    inline sint16 card() const { return set.size(); }

    void print(IndentedOStream &os = ind_cout) const {
        os << "There " << (card() > 1 ? "are " : "is ") << card() << " element"
           << (card() > 1 ? "s " : " ") << "in the super summit set."
           << EndLine(2);

        os << "─────" << EndLine() << " Set " << EndLine() << "─────";

        os.Indent(4);

        os << EndLine(1);

        sint16 indent = (int(std::to_string(card() - 1).length()) + 1) / 4 + 1;
        sint16 count = 0;
        for (ConstIterator it = begin(); it != end(); it++) {
            os << count << ":";
            for (sint16 _ = 0;
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

template <class F>
SuperSummitSet<BraidTemplate<F>> super_summit_set(const BraidTemplate<F> &b) {
    std::list<BraidTemplate<F>> queue, queue_rcf;
    SuperSummitSet<BraidTemplate<F>> sss;

    BraidTemplate<F> b2 = send_to_super_summit(b);
    BraidTemplate<F> b2_rcf = b2;
    b2_rcf.MakeRCFFromLCF();

    queue.push_back(b2);
    queue_rcf.push_back(b2_rcf);

    sss.insert(b2);

    while (!queue.empty()) {
        std::vector<F> min = min_super_summit(queue.front(), queue_rcf.front());

        for (typename std::vector<F>::iterator itf = min.begin();
             itf != min.end(); itf++) {
            b2 = queue.front();
            b2.Conjugate(*itf);

            if (!sss.mem(b2)) {
                b2_rcf = queue_rcf.front();
                b2_rcf.ConjugateRCF(*itf);

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

template <class F>
inline bool are_conjugate(const BraidTemplate<F> &u,
                          const BraidTemplate<F> &v) {
    std::unordered_set<BraidTemplate<F>> u_sss = super_summit_set(u);
    return u_sss.mem(send_to_super_summit(v));
}

} // namespace cgarside::super_summit

#endif
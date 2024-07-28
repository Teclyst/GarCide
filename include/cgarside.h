#ifndef CGARSIDE
#define CGARSIDE

#include "utility.h"
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Generic factors and braids representation for Garside groups.
 *
 * In this namespace, we define (templated) classes that represent elements for
 * generic Garside groups.
 */
namespace cgarside {

template <class U> class FactorTemplate {

  public:
    using Parameter = typename U::Parameter;

  private:
    // underlying is the data structure that actually represents the factor
    // (e.g., a permutation table for a braid canonical factor).
    U underlying;

  public:
    // FactorTemplate(under) initializes a new factor, with underlying element
    // under.
    FactorTemplate(const U &under) : underlying(under) {}

    FactorTemplate(Parameter parameter) : underlying(parameter) {}

    inline U get_underlying() const { return underlying; }

    inline static Parameter parameter_of_string(const std::string &str) {
        return U::parameter_of_string(str);
    }

    inline Parameter get_parameter() const {
        return underlying.get_parameter();
    }

    inline sint16 lattice_height() const {
        return underlying.lattice_height();
    };

    // a.of_string sets a to the factor specified by str.
    void of_string(const std::string &str, size_t &pos) {
        underlying.of_string(str, pos);
    }

    // a.debug(os) prints a's internal representation to os.
    void debug(IndentedOStream &os = ind_cout) const {
        os << "{   Underlying:";
        os.Indent(8);
        os << EndLine();
        underlying.debug(os);
        os.Indent(-8);
        os << EndLine() << "}";
    }

    // a.print(os) prints a to os.
    void print(IndentedOStream &os = ind_cout) const { underlying.print(os); }

    // a.identity sets a to identity.
    void identity() { underlying.identity(); }

    // a.delta() sets a to delta.
    void delta() { underlying.delta(); }

    // a.compare(b) returns true if a and b are equal, false otherwise.
    bool compare(const FactorTemplate &b) const {
        return underlying.compare(b.underlying);
    }

    // a == b returns true if a and b are equal, false otherwise.
    // Syntactic sugar for a.compare(b).
    bool operator==(const FactorTemplate &b) const { return compare(b); }

    // a != b returns true if a and b are not equal, false otherwise.
    bool operator!=(const FactorTemplate &b) const { return !compare(b); }

    // a.is_delta() returns whether a == e.
    bool is_identity() const {
        FactorTemplate e = FactorTemplate(*this);
        e.identity();
        return compare(e);
    }

    // a.is_delta() returns whether a = delta.
    bool is_delta() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.delta();
        return compare(delta);
    }

    // a.left_complement(b) returns (assuming that a right-divides b) the left
    // complement of a under b, ba^{-1}.
    FactorTemplate left_complement(const FactorTemplate &b) const {
        return FactorTemplate(underlying.left_complement(b.underlying));
    }

    // a.left_complement() return a's left complement.
    FactorTemplate left_complement() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.delta();
        return left_complement(delta);
    }

    // a.right_complement(b) returns (assuming that a left-divides b) the right
    // complement of a under b, a^{-1}b.
    FactorTemplate right_complement(const FactorTemplate &b) const {
        return FactorTemplate(underlying.right_complement(b.underlying));
    }

    // a.right_complement() return a's right complement.
    FactorTemplate right_complement() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.delta();
        return right_complement(delta);
    }

    // ~a return a's right complement.
    // Syntactic sugar for a.right_complement().
    FactorTemplate operator~() const { return right_complement(); }

    // Syntactic sugar for b.right_complement(a).
    FactorTemplate operator/(const FactorTemplate &b) const {
        return b.right_complement(*this);
    }

    void delta_conjugate_mut(sint16 k) { underlying.delta_conjugate_mut(k); }

    // a.delta_conjugate(k) returns a, conjugated by Delta ^ k.
    // Makes 2 |k| complement calculations.
    FactorTemplate delta_conjugate(sint16 k) const {
        FactorTemplate conjugate = *this;
        conjugate.delta_conjugate_mut(k);
        return conjugate;
    }

    // a.delta_conjugate() returns a conjugated by Delta.
    FactorTemplate delta_conjugate() const { return delta_conjugate(1); }

    // a.left_meet(b) returns the left meet of a and b.
    FactorTemplate left_meet(const FactorTemplate &b) const {
        return FactorTemplate(underlying.left_meet(b.underlying));
    }

    // a ^ b returns the left meet of a and b.
    // Syntactic sugar for a.left_meet(b).
    FactorTemplate operator^(const FactorTemplate &b) const {
        return left_meet(b);
    }

    // a.right_meet(b) returns the right meet of a and b.
    FactorTemplate right_meet(const FactorTemplate &b) const {
        return FactorTemplate(underlying.right_meet(b.underlying));
    }

    // a.left_join(b) returns the left join of a and b.
    FactorTemplate left_join(const FactorTemplate &b) const {
        return right_complement()
            .right_meet(b.right_complement())
            .left_complement();
    }

    // a.right_join(b) returns the right join of a and b.
    FactorTemplate right_join(const FactorTemplate &b) const {
        return left_complement().left_meet(b.left_complement()).right_complement();
    }

    // a.is_left_weighted(b) returns true if a | b is left weighted, or false
    // otherwise.
    bool is_left_weighted(const FactorTemplate &b) const {
        return right_complement().left_meet(b).is_identity();
    }

    // a.IsRightWeighted(b) returns true if a | b is right weighted, or false
    // otherwise.
    bool IsRightWeighted(const FactorTemplate &b) const {
        return left_meet(b.left_complement()).is_identity();
    }

    // a.product(b) returns the product of two factors, under the assumption
    // that it lies below Delta.
    FactorTemplate product(const FactorTemplate &b) const {
        return FactorTemplate(underlying.product(b.underlying));
    }

    void right_multiply(const FactorTemplate &b) { *this = *this * b; }

    std::size_t hash() const { return underlying.hash(); }

    // a * b is the product of a and b, under the assumption that it lies below
    // Delta. Syntactic sugar for a.product(b).
    FactorTemplate operator*(const FactorTemplate &b) const {
        return product(b);
    }

    // a.randomize() sets a to a random factor.
    void randomize() {

#ifdef RANDOMIZE_ON_ATOMS

        std::vector<FactorTemplate> atoms = atoms();

        *this = atoms[rand() % int(atoms.size())];

#else

        underlying.randomize();

#endif
    }

    // a.atoms() returns the list of the atoms.
    std::vector<FactorTemplate> atoms() const {
        std::vector<U> atoms = underlying.atoms();
        typename std::vector<U>::iterator atoms_it;
        std::vector<FactorTemplate> factor_atoms;
        for (auto const &atoms_it : atoms) {
            factor_atoms.push_back(atoms_it);
        }
        return factor_atoms;
    }
};

// make_left_weighted(u, v) computes the left-weighted decomposition u' | v' =
// u | v, and sets u = u' and v = v'. It then returns true if something was
// done (so that it may be used with `apply_binfun`). SHOULD NEVER BE CALLED
// UPON u, v IF
// `&u == &v`!
template <class F> bool make_left_weighted(F &u, F &v) {
    F t = (~u) ^ v;
    if (t.is_identity()) {
        return false;
    } else {
        v = v / t;
        u = u * t;
        return true;
    }
}

// make_right_weighted(u, v) computes the right-weighted decomposition u' | v'
// = u | v, and sets u = u' and v = v'. It then returns true if something
// was done (so that it may be used with `apply_binfun`). SHOULD NEVER BE
// CALLED UPON u, v IF `&u == &v`!
template <class F> bool make_right_weighted(F &u, F &v) {
    F t = u.right_meet(v.left_complement());
    if (t.is_identity()) {
        return false;
    } else {
        v = t * v;
        u = t.left_complement(u);
        return true;
    }
}

// Overloading << for factor classes.
template <class U>
IndentedOStream &operator<<(IndentedOStream &os, const FactorTemplate<U> &f) {
    f.print(os);
    return os;
}

/**
 * @brief A class representing braids for generic Garside groups.
 *
 * This class represents elements of a Garside groups, whose canonicals
 * factors are represented by the class `F`.
 *
 * Braids are maintained in LCF most of the time. It should be assumed that
 * a method that does not have "RCF" in its name expects all of its
 * arguments to be in LCF.
 *
 * @tparam F A class that represents factors. It does not have to be a
 * `FactorTemplate<U>` for some class `U`, but has to implement all methods of
 * `FactorTemplate` classes.
 */
template <class F> class BraidTemplate {

  public:
    /**
     * @brief A type for parameters.
     *
     * We often implement Garside structures for a series of groups, that
     * are distinguished one from another by some parameter (often, but not
     * always, an integer).
     *
     */
    using Parameter = typename F::Parameter;

    /**
     * @brief A (group) parameter.
     *
     * We often implement Garside structures for a series of groups, that
     * are distinguished one from another by some parameter (often, but not
     * always, an integer).
     *
     */
    Parameter Parameter;

    // `Delta` is the number of Deltas on the left end of the word.

    /**
     * @brief Infimum.
     *
     * The power of delta at the left end of the word (in LCF, otherwise at
     * the right end in RCF).
     *
     */
    sint32 Delta;

    /**
     * @brief The braid's canonical factors.
     *
     * A list of the braid's canonical factors, from left to right. It is
     * left weighted when in LCF, and right weighted when in RCF.
     */
    std::list<F> FactorList;

  public:
    using FactorItr = typename std::list<F>::iterator;
    using RevFactorItr = typename std::list<F>::reverse_iterator;
    using ConstFactorItr = typename std::list<F>::const_iterator;
    using ConstRevFactorItr = typename std::list<F>::const_reverse_iterator;

    inline FactorItr begin() const {
        return FactorList.begin();
    }

    inline RevFactorItr rbegin() const {
        return FactorList.rbegin();
    }

    inline ConstFactorItr cbegin() const {
        return FactorList.begin();
    }

    inline ConstRevFactorItr crbegin() const {
        return FactorList.rbegin();
    }

    inline FactorItr end() const {
        return FactorList.end();
    }

    inline RevFactorItr rend() const {
        return FactorList.rend();
    }

    inline ConstFactorItr cend() const {
        return FactorList.end();
    }

    inline ConstRevFactorItr crend() const {
        return FactorList.rend();
    }

  public:
    /**
     * @brief Construct a new BraidTemplate, with a group parameter.
     *
     * Construct a new BraidTemplate, with `Parameter = parameter`. It is
     * initialized as the identity braid.
     *
     * @param parameter Group parameter.
     */
    BraidTemplate(Parameter parameter)
        : Parameter(parameter), Delta(0), FactorList() {}

    /**
     * @brief Construct a new BraidTemplate, from a factor.
     *
     * Construct a new BraidTemplate whose only factor is `f`.
     *
     * @param f FactorTemplate to be converted to a braid.
     */
    BraidTemplate(const F &f)
        : Parameter(f.get_parameter()), Delta(0), FactorList() {
        if (f.is_delta()) {
            Delta = 1;
        } else if (!f.is_identity()) {
            FactorList.push_back(f);
        }
    }

    inline static Parameter parameter_of_string(const std::string &str) {
        return F::parameter_of_string(str);
    }

    Parameter get_parameter() const { return Parameter; }

    /**
     * @brief Prints `*this` to `os`.
     *
     * Prints `*this` to `os`, in a format that is compatible with
     * `of_string` (assuming that the same holds for `F`).
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os = ind_cout) const {
        if (Delta != 0 && Delta != 1) {
            os << "D ^ " << Delta << (CanonicalLength() > 0 ? " . " : "");
        } else if (Delta == 1) {
            os << "D" << (CanonicalLength() > 0 ? " . " : "");
        }
        ConstFactorItr it_end = cend();
        it_end--;
        for (ConstFactorItr it = cbegin(); it != cend();
             it++) {
            (*it).print(os);
            if (it != it_end) {
                os << " . ";
            }
        }
    }

    /**
     * @brief Prints `*this` (in RCF) to `os`.
     *
     * Prints `*this` to `os`, in a format that is compatible with
     * `of_string` (assuming that the same holds for `F`).
     *
     * @param os The output stream it prints to.
     */
    void print_rcf(IndentedOStream &os = ind_cout) const {
        ConstFactorItr it_end = FactorList.end();
        it_end--;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            (*it).print(os);
            if (it != it_end) {
                os << " . ";
            }
        }
        if (Delta != 0 && Delta != 1) {
            os << (CanonicalLength() > 0 ? " . " : "") << "D ^ " << Delta;
        } else if (Delta == 1) {
            os << (CanonicalLength() > 0 ? " . " : "") << "D";
        }
    }

    // `w.identity()` sets w to the empty word.
    inline void identity() {
        Delta = 0;
        FactorList.clear();
    }

    // `u.CanonicalLength` returns u's canonical length.
    inline sint16 CanonicalLength() const { return FactorList.size(); }

    inline sint16 Inf() const { return Delta; }

    inline sint16 Sup() const { return Inf() + CanonicalLength(); }

    // `u.compare(v)` returns whether u and v have the same internal
    // representation.
    inline bool compare(const BraidTemplate &v) const {
        return (Delta == v.Delta && FactorList == v.FactorList);
    }

    // `u == v` returns whether u and v have the same internal
    // representation. Syntactic sugar for `u.compare(v)`.
    bool operator==(const BraidTemplate &v) const { return compare(v); }

    // `u != v` returns whether u and v do not have the same internal
    // representation.
    bool operator!=(const BraidTemplate &v) const { return !compare(v); }

    // `u.is_identity` returns whether u represents the identity element.
    bool is_identity() const { return Delta == 0 && FactorList.empty(); }

    // `u.Inverse()` returns the inverse of u.
    //  See the ElRifai and Morton 1994 article for correction.
    BraidTemplate Inverse() const {
        BraidTemplate b(get_parameter());
        b.Delta = -Delta;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            // Rewrite a_1 ... a_k (f)^(- 1) Delta^r as
            // a_1 ... a_k Delta ^ (r - 1) (Delta^(- r) d_L(f) Delta^r).
            b.FactorList.push_front(
                (*it).left_complement().delta_conjugate(b.Delta));
            --b.Delta;
        }
        return b;
    }

    // `u.Inverse()` returns the inverse of u.
    //  See the ElRifai and Morton 1994 article for correction.
    BraidTemplate InverseRCF() const {
        BraidTemplate b(get_parameter());
        b.Delta = -Delta;
        for (ConstRevFactorItr revit = FactorList.rbegin();
             revit != FactorList.rend(); revit++) {
            // Rewrite Delta^r (f)^(- 1) a_1 ... a_k as
            // (Delta^r d_R(f) Delta^(- r)) Delta ^ (r - 1) a_1 ... a_k.
            b.FactorList.push_back(
                (*revit).right_complement().delta_conjugate(-b.Delta));
            --b.Delta;
        }
        return b;
    }

    // `!u` returns the inverse of u.
    // Syntactic sugar for `u.Inverse()`.
    BraidTemplate operator!() const { return Inverse(); }

    // Clean gets rid of (factor) Deltas at the start, and identity elements
    // at the end pf `FactorList`.
    void Clean() {
        FactorItr it = FactorList.begin();
        while (it != FactorList.end() && (*it).is_delta()) {
            ++it;
            ++Delta;
        }
        FactorList.erase(FactorList.begin(), it);
        RevFactorItr revit = FactorList.rbegin();
        while (revit != FactorList.rend() && (*revit).is_identity()) {
            ++revit;
        }
        FactorList.erase(revit.base(), FactorList.end());
    }

    void CleanRCF() {
        FactorItr it = FactorList.begin();
        while (it != FactorList.end() && (*it).is_identity()) {
            ++it;
        }
        FactorList.erase(FactorList.begin(), it);
        RevFactorItr revit = FactorList.rbegin();
        while (revit != FactorList.rend() && (*revit).is_delta()) {
            ++revit;
            ++Delta;
        }
        FactorList.erase(revit.base(), FactorList.end());
    }

    // `u.LeftProduct(f)` assigns fu to u.
    void LeftProduct(const F &f) {
        FactorList.push_front(f.delta_conjugate(Delta));
        apply_binfun(FactorList.begin(), FactorList.end(), make_left_weighted<F>);
        Clean();
    }

    // `u.right_multiply(f)` assigns uf to u.
    void right_multiply(const F &f) {
        FactorList.push_back(f);
        reverse_apply_binfun(FactorList.begin(), FactorList.end(),
                             make_left_weighted<F>);
        Clean();
    }

    // `u.LeftProduct(v)` assigns v u to u.
    void LeftProduct(const BraidTemplate &v) {
        for (ConstRevFactorItr it = v.FactorList.rbegin();
             it != v.FactorList.rend(); it++) {
            LeftProduct(*it);
        }
        Delta += v.Delta;
    }

    // `u.right_multiply(v)` assigns u v to u.
    // v's factors move directly to u - be careful.
    void right_multiply(const BraidTemplate &v) {
        for (FactorItr it = FactorList.begin(); it != FactorList.end(); it++) {
            (*it) = (*it).delta_conjugate(v.Delta);
        }
        Delta += v.Delta;
        for (ConstFactorItr it = v.FactorList.begin(); it != v.FactorList.end();
             it++) {
            right_multiply((*it));
        }
    }

    // `u.LeftDivide(v)` assigns v ^ (- 1) u to u.
    inline void LeftDivide(const BraidTemplate &v) { LeftProduct(!v); }

    // `u.RightDivide(v)` assigns u v ^ (- 1) to u.
    inline void RightDivide(const BraidTemplate &v) { right_multiply(!v); }

    // `u.LeftDivide(f)` assigns f ^ (- 1) u to u.
    inline void LeftDivide(const F &f) { LeftProduct(!BraidTemplate(f)); }

    // `u.RightDivide(v)` assigns u v ^ (- 1) to u.
    inline void RightDivide(const F &f) { right_multiply(!BraidTemplate(f)); }

    // `u.LeftProductRCF(f)` assigns fu to u.
    void LeftProductRCF(const F &f) {
        FactorList.push_front(f);
        apply_binfun(FactorList.begin(), FactorList.end(),
                     make_right_weighted<F>);
        CleanRCF();
    }

    // `u.right_multiply(f)` assigns uf to u.
    void RightProductRCF(const F &f) {
        FactorList.push_back(f.delta_conjugate(-Delta));
        reverse_apply_binfun(FactorList.begin(), FactorList.end(),
                             make_right_weighted<F>);
        CleanRCF();
    }

    // `u.LeftProduct(v)` assigns v u to u.
    void LeftProductRCF(const BraidTemplate &v) {
        for (RevFactorItr it = FactorList.rbegin(); it != FactorList.rend();
             it++) {
            (*it) = (*it).delta_conjugate(-v.Delta);
        }
        Delta += v.Delta;
        for (ConstRevFactorItr it = v.FactorList.rbegin();
             it != v.FactorList.rend(); it++) {
            LeftProductRCF(*it);
        }
    }

    // `u.right_multiply(v)` assigns u v to u.
    void RightProductRCF(const BraidTemplate &v) {
        for (ConstFactorItr it = v.FactorList.begin(); it != v.FactorList.end();
             it++) {
            RightProductRCF(*it);
        }
        Delta += v.Delta;
    }

    BraidTemplate left_meet(const BraidTemplate &v) const {
        sint16 shift = 0;
        BraidTemplate b = BraidTemplate(get_parameter());
        F f1 = F(get_parameter()), f2 = F(get_parameter()),
          f = F(get_parameter());
        f.Delta();

        BraidTemplate b1 = *this, b2 = v;

        shift -= b1.Delta;
        b2.Delta -= b1.Delta;
        b1.Delta = 0;

        if (b2.Delta < 0) {
            shift -= b2.Delta;
            b1.Delta -= b2.Delta;
            b2.Delta = 0;
        }

        while (!f.is_identity()) {
            if (b1.Delta > 0) {
                f1.Delta();
            } else if (b1.CanonicalLength() == 0) {
                f1.identity();
            } else {
                f1 = b1.FactorList.front();
            }

            if (b2.Delta > 0) {
                f2.Delta();
            } else if (b2.CanonicalLength() == 0) {
                f2.identity();
            } else {
                f2 = b2.FactorList.front();
            }

            f = f1 ^ f2;

            b.right_multiply(f);
            b1.LeftDivide(f);
            b2.LeftDivide(f);
        }

        b.Delta -= shift;
        return b;
    }

    inline BraidTemplate left_meet(const F &f) {
        return left_meet(BraidTemplate(f));
    }

    inline BraidTemplate operator^(const BraidTemplate &v) {
        return left_meet(v);
    }

    inline BraidTemplate operator^(const F &f) { return left_meet(f); }

    BraidTemplate left_join(const BraidTemplate &v) const {
        sint16 shift = 0;
        BraidTemplate b = BraidTemplate(get_parameter());
        F f2 = F(get_parameter()), f = F(get_parameter());
        f.Delta();

        BraidTemplate b1 = *this, b2 = v;

        shift -= b1.Delta;
        b2.Delta -= b1.Delta;
        b1.Delta = 0;

        if (b2.Delta < 0) {
            shift -= b2.Delta;
            b1.Delta -= b2.Delta;
            b2.Delta = 0;
        }

        b = b1;

        while (!b2.is_identity()) {
            if (b2.Delta > 0) {
                f2.Delta();
            } else if (b2.CanonicalLength() == 0) {
                f2.identity();
            } else {
                f2 = b2.FactorList.front();
            }

            f = b1.Remainder(f2);

            b.right_multiply(f);
            b1.right_multiply(f);
            b1.LeftDivide(f2);
            b2.LeftDivide(f2);
        }

        b.Delta -= shift;
        return b;
    }

    inline BraidTemplate left_join(const F &f) const {
        return left_join(BraidTemplate(f));
    }

    inline BraidTemplate right_meet(const BraidTemplate &v) const {
        return !((!(*this)).left_join(!v));
    }

    inline BraidTemplate right_meet(const F &f) const {
        return !((!(*this)).left_join(!BraidTemplate(f)));
    }

    inline BraidTemplate right_join(const BraidTemplate &v) const {
        return !((!(*this)).left_meet(!v));
    }

    inline BraidTemplate right_join(const F &f) const {
        return !((!(*this)).left_meet(!BraidTemplate(f)));
    }

    // `u.LeftDivide(v)` assigns v ^ (- 1) u to u.
    inline void LeftDivideRCF(const BraidTemplate &v) {
        LeftProductRCF(v.InverseRCF());
    }

    // `u.RightDivide(v)` assigns u v ^ (- 1) to u.
    inline void RightDivideRCF(const BraidTemplate &v) {
        RightProductRCF(v.InverseRCF());
    }

    // `u.LeftDivide(f)` assigns f ^ (- 1) u to u.
    inline void LeftDivideRCF(const F &f) {
        LeftProductRCF(BraidTemplate(f).InverseRCF());
    }

    // `u.RightDivide(v)` assigns u v ^ (- 1) to u.
    inline void RightDivideRCF(const F &f) {
        RightProductRCF(BraidTemplate(f).InverseRCF());
    }

    inline void Conjugate(const F &f) {
        LeftDivide(f);
        right_multiply(f);
    }

    inline void Conjugate(const BraidTemplate &v) {
        LeftDivide(v);
        right_multiply(v);
    }

    inline void ConjugateRCF(const F &f) {
        LeftDivideRCF(f);
        RightProductRCF(f);
    }

    inline void ConjugateRCF(const BraidTemplate &v) {
        LeftDivideRCF(v);
        RightProductRCF(v);
    }

    // `u.Initial()` returns the initial factor of u, that is, if u = Delta
    // ^ r u_1 ... u_k, Delta ^ r u_1 Delta ^ (- r). If u has canonical
    // length zero, returns the identity factor instead.
    inline F Initial() const {
        if (CanonicalLength() == 0) {
            F id = F(get_parameter());
            id.identity();
            return id;
        } else {
            return FactorList.front().delta_conjugate(-Delta);
        }
    }

    inline F First() const {
        if (CanonicalLength() == 0) {
            F id = F(get_parameter());
            id.identity();
            return id;
        } else {
            return FactorList.front();
        }
    }

    // `u.Final()` returns the final factor of u, that is, if u = Delta ^ r
    // u_1
    // ... u_k, u_k. If u has canonical length zero, returns the identity
    // factor instead.
    inline F Final() const {
        if (CanonicalLength() == 0) {
            F id = F(get_parameter());
            id.identity();
            return id;
        } else {
            return FactorList.back();
        }
    }

    // `u.PreferredPrefix()` returns the preferred prefix of u, that is, if
    // u = Delta ^ r u_1 ... u_k, p(u) = d_R(u_k) ^_L Delta ^ r u_1 Delta ^
    // (- r) If u has canonical length zero, returns the identity factor
    // instead.
    inline F PreferredPrefix() const { return Initial() ^ ~Final(); }

    F PreferredSuffixRCF() const {
        if (CanonicalLength() == 0) {
            F id = F(get_parameter());
            id.identity();
            return id;
        } else {
            return FactorList.back().delta_conjugate(Delta).right_meet(
                FactorList.front().left_complement());
        }
    };

    F PreferredSuffix() const {
        BraidTemplate right = BraidTemplate(*this);
        right.MakeRCFFromLCF();
        return right.PreferredSuffixRCF();
    };

    // `u.Cycling()` cycles u: if u = Delta ^ r u_1 ... u_k, then after
    // applying cycling u will contain (the LNF of) Delta ^ r u_2 ... u_k
    // (Delta ^ r u_1 Delta ^ (-r)).
    inline void Cycling() {
        if (CanonicalLength() == 0) {
            return;
        }
        F i = Initial();
        FactorList.pop_front();
        right_multiply(i);
    }

    // `u.Decycling()` decycles u: if u = Delta ^ r u_1 ... u_k, then after
    // applying cycling u will contain (the LNF of) Delta ^ r u_2 ... u_k
    // (Delta ^ r u_1 Delta ^ (-r)).
    inline void Decycling() {
        if (CanonicalLength() == 0) {
            return;
        }
        F f = Final();
        FactorList.pop_back();
        LeftProduct(f);
    }

    // `u.Sliding()` cyclically slides u: if u = Delta ^ r u_1 ... u_k, and
    // Delta ^ r u_1 Delta ^ (-r) = p(u) u'_1, then after applying cycling u
    // will contain (the LNF of) Delta ^ r (Delta ^ r u'_1 Delta ^ (-r)) u_2
    // ... u_k p(u).
    inline void Sliding() {
        if (CanonicalLength() == 0) {
            return;
        }
        Conjugate(PreferredPrefix());
    }

    // `u.product(v)` returns uv.
    BraidTemplate product(const BraidTemplate &v) const {
        BraidTemplate w(*this);
        w.right_multiply(v);
        return w;
    }

    // `u.Power(k)` returns u raised to the power k.
    // Uses a fast exponentiation algorithm; the number of multiplications
    // is logarithmic in k.
    BraidTemplate Power(const sint16 k) const {
        if (k == 0) {
            return BraidTemplate(get_parameter());
        } else if (k % 2 == 0) {
            BraidTemplate root = Power(k / 2);
            return root * root;
        } else if (k > 0) {
            BraidTemplate root = Power(k / 2);
            return *this * root * root;
        } else {
            BraidTemplate root = Power(k / 2);
            return !*this * root * root;
        }
    }

    // `u * v` returns uv.
    // Syntactic sugar for `u.product(v)`.
    BraidTemplate operator*(const BraidTemplate &v) const { return product(v); }

    // `u.Normalize()` turns u into LCF.
    inline void Normalize() {
        bubble_sort(FactorList.begin(), FactorList.end(), make_left_weighted<F>);
        Clean();
    }

    // `u.LCFToRCF()` turns u, assumed to be in LCF (so that we may avoid
    // having to clean up), into RCF. The result (r, [u_1, ..., u_k])
    // represents u_1
    // ... u_k Delta ^ r.
    inline void MakeRCFFromLCF() {
        for (FactorItr it = FactorList.begin(); it != FactorList.end(); ++it) {
            *it = (*it).delta_conjugate(-Delta);
        };
        bubble_sort(FactorList.begin(), FactorList.end(), make_right_weighted<F>);
    }

    // `u.LCFToRCF()` turns u, assumed to be in LCF, into RCF (so that we
    // may avoid having to clean up). The result (r, [u_1, ..., u_k])
    // represents u_1
    // ... u_k Delta ^ r.
    inline void MakeLCFFromRCF() {
        for (FactorItr it = FactorList.begin(); it != FactorList.end(); ++it) {
            *it = (*it).delta_conjugate(Delta);
        };
        bubble_sort(FactorList.begin(), FactorList.end(), make_left_weighted<F>);
    }

    // `b.Remainder(f)` computes, if b is positive, the simple factor s such
    // that bs is the left lcm of b and b and f.
    F Remainder(const F &f) const {
        F fi = f;
        if (Delta != 0) {
            fi.identity();
        } else {
            for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
                 it++) {
                fi = (*it).left_join(fi) / *it;
            }
        }
        return fi;
    }

    sint16 Rigidity() const {
        BraidTemplate b2 = *this;
        sint16 rigidity = 0;

        if (CanonicalLength() == 0)
            return rigidity;

        b2.right_multiply(b2.Initial());

        ConstFactorItr it2 = b2.FactorList.begin();

        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++, it2++, rigidity++) {
            if (*it != *it2) {
                break;
            }
        }

        return rigidity;
    }

    // Randomizes the braid, setting it at a given length, using parameters
    // extracted from FactorTemplate f. The result isn't in LCF.
    void randomize(sint16 canonical_length) {
        FactorList.clear(); // Potential memory leak? To be checked.
        Delta = 0;
        for (sint16 i = 0; i < canonical_length; i++) {
            F f = F(get_parameter());
            f.randomize();
            FactorList.push_back(f);
        }
    }

    void debug(IndentedOStream &os = ind_cout) const {
        os << "{   ";
        os.Indent(4);
        os << "Parameter:";
        os.Indent(4);
        os << EndLine();
        os << Parameter;
        os.Indent(-4);
        os << EndLine();
        os << "Delta:";
        os.Indent(4);
        os << EndLine();
        os << Delta;
        os.Indent(-4);
        os << EndLine();
        os << "FactorList:";
        os.Indent(4);
        os << EndLine();
        os << "[   ";
        os.Indent(4);
        ConstFactorItr it_end = FactorList.end();
        it_end--;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            (*it).debug(os);
            if (it != it_end) {
                os << "," << EndLine();
            }
        }
        os.Indent(-4);
        os << EndLine() << "]";
        os.Indent(-8);
        os << EndLine() << "}";
    }

    std::size_t hash() const {
        std::size_t h = Delta;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            h = h * 31 + (*it).hash();
        }
        return h;
    }

    /**
     * @brief Conversion from string.
     *
     * Reads the string `str` and sets `this` to the corresponding braid, in
     * LCF.
     *
     * Letting `L` be the language of (a generating subset of) the factors,
     * `W =
     * (\s | \t)*` be the language of whitespaces, and `Z = -? ([1 - 9] [0 -
     * 9]* | 0)` be the language of integers, accepted strings are those
     * represented by regular expression `(W | .)* (((D | L) (W ^ W Z)?) (W
     * | .)*)*`
     *
     * @param str The string to convert from.
     * @exception `InvalidStringError`: Thrown when it isn't possible to
     * extract a factor from `str`, or when the factor we tried to extract
     * does not exist (e.g. `4` isn't a legal factor for artin braids on 4
     * strands).
     */
    void of_string(const std::string str) {
        size_t pos = 0;

        // We use double backslashes, as we want to obtain escape sequences.
        // Otherwise, for instance, `\t` would be recognized as a tab by
        // C++. For another example, `\^` would be recognized as an illegal
        // escape sequence, and compilation would fail. If we wanted to
        // recognize a single backslash `\`, we would have to escape it as
        // `\\\\`. In C++ `std::regex`, `\s` is a whitespace, `\t` a tab,
        // `\.` a dot and `\^` a chevron (`.` and `^` are special characters
        // for `std::regex`).
        std::regex ignore{"[\\s\\.\\t]*"};
        std::regex power{"[\\s\\t]*\\^[\\s\\t]*(" + number_regex + ")"};
        std::regex inverse{"![\\s\\t]*"};

        BraidTemplate b = *this;
        b.identity();

        std::smatch match;

        F fact(get_parameter());

        std::regex_search(str.begin() + pos, str.end(), match, ignore,
                          std::regex_constants::match_continuous);
        pos += match[0].length();

        while (pos != str.length()) {
            sint16 pow;

            fact.of_string(str, pos);

            if (std::regex_search(str.begin() + pos, str.end(), match, power,
                                  std::regex_constants::match_continuous)) {
                pos += match[0].length();
                pow = std::stoi(match[1]);
            } else {
                pow = 1;
            }
            if (pow >= 0) {
                for (sint16 _ = 0; _ < pow; _++) {
                    b.right_multiply(fact);
                }
            } else {
                pow = -pow;
                for (sint16 _ = 0; _ < pow; _++) {
                    b.RightDivide(fact);
                }
            }
            std::regex_search(str.begin() + pos, str.end(), match, ignore,
                              std::regex_constants::match_continuous);
            pos += match[0].length();
        }
        *this = b;
    }
};

// Overloading << for braid classes.
template <class F>
IndentedOStream &operator<<(IndentedOStream &os, const BraidTemplate<F> &b) {
    b.print(os);
    return os;
}

} // namespace cgarside

template <class U> struct std::hash<cgarside::FactorTemplate<U>> {
    std::size_t
    operator()(cgarside::FactorTemplate<U> const &f) const noexcept {
        return f.hash();
    }
};

template <class F> struct std::hash<cgarside::BraidTemplate<F>> {
    std::size_t operator()(cgarside::BraidTemplate<F> const &u) const noexcept {
        return u.hash();
    }
};

#endif
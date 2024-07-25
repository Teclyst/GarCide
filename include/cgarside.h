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
    typedef typename U::ParameterType ParameterType;

  private:
    // Underlying is the data structure that actually represents the factor
    // (e.g., a permutation table for a braid canonical factor).
    U Underlying;

  public:
    // FactorTemplate(under) initializes a new factor, with underlying element under.
    FactorTemplate(const U &under) : Underlying(under) {}

    FactorTemplate(ParameterType parameter) : Underlying(parameter) {}

    U GetUnderlying() const { return Underlying; }

    inline static ParameterType parameter_of_string(const std::string &str) {
        return U::parameter_of_string(str);
    }

    inline ParameterType GetParameter() const {
        return Underlying.GetParameter();
    }

    inline sint16 LatticeHeight() const { return Underlying.LatticeHeight(); };

    // a.OfString sets a to the factor specified by str.
    void OfString(const std::string &str, size_t &pos) {
        Underlying.OfString(str, pos);
    }

    // a.Debug(os) prints a's internal representation to os.
    void Debug(IndentedOStream &os = ind_cout) const {
        os << "{   Underlying:";
        os.Indent(8);
        os << EndLine();
        Underlying.Debug(os);
        os.Indent(-8);
        os << EndLine() << "}";
    }

    // a.Print(os) prints a to os.
    void Print(IndentedOStream &os = ind_cout) const { Underlying.Print(os); }

    // a.Identity sets a to Identity.
    void Identity() { Underlying.Identity(); }

    // a.Delta() sets a to Delta.
    void Delta() { Underlying.Delta(); }

    // a.Compare(b) returns true if a and b are equal, false otherwise.
    bool Compare(const FactorTemplate &b) const {
        return Underlying.Compare(b.Underlying);
    }

    // a == b returns true if a and b are equal, false otherwise.
    // Syntactic sugar for a.Compare(b).
    bool operator==(const FactorTemplate &b) const { return Compare(b); }

    // a != b returns true if a and b are not equal, false otherwise.
    bool operator!=(const FactorTemplate &b) const { return !Compare(b); }

    // a.IsDelta() returns whether a == e.
    bool IsIdentity() const {
        FactorTemplate e = FactorTemplate(*this);
        e.Identity();
        return Compare(e);
    }

    // a.IsDelta() returns whether a = Delta.
    bool IsDelta() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.Delta();
        return Compare(delta);
    }

    // a.LeftComplement(b) returns (assuming that a right-divides b) the left
    // complement of a under b, ba^{-1}.
    FactorTemplate LeftComplement(const FactorTemplate &b) const {
        return FactorTemplate(Underlying.LeftComplement(b.Underlying));
    }

    // a.LeftComplement() return a's left complement.
    FactorTemplate LeftComplement() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.Delta();
        return LeftComplement(delta);
    }

    // a.RightComplement(b) returns (assuming that a left-divides b) the right
    // complement of a under b, a^{-1}b.
    FactorTemplate RightComplement(const FactorTemplate &b) const {
        return FactorTemplate(Underlying.RightComplement(b.Underlying));
    }

    // a.RightComplement() return a's right complement.
    FactorTemplate RightComplement() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.Delta();
        return RightComplement(delta);
    }

    // ~a return a's right complement.
    // Syntactic sugar for a.RightComplement().
    FactorTemplate operator~() const { return RightComplement(); }

    // Syntactic sugar for b.RightComplement(a).
    FactorTemplate operator/(const FactorTemplate &b) const { return b.RightComplement(*this); }

    void DeltaConjugateMut(sint16 k) { Underlying.DeltaConjugate(k); }

    // a.DeltaConjugate(k) returns a, conjugated by Delta ^ k.
    // Makes 2 |k| complement calculations.
    FactorTemplate DeltaConjugate(sint16 k) const {
        FactorTemplate conjugate = *this;
        conjugate.DeltaConjugateMut(k);
        return conjugate;
    }

    // a.DeltaConjugate() returns a conjugated by Delta.
    FactorTemplate DeltaConjugate() const { return DeltaConjugate(1); }

    // a.LeftMeet(b) returns the left meet of a and b.
    FactorTemplate LeftMeet(const FactorTemplate &b) const {
        return FactorTemplate(Underlying.LeftMeet(b.Underlying));
    }

    // a ^ b returns the left meet of a and b.
    // Syntactic sugar for a.LeftMeet(b).
    FactorTemplate operator^(const FactorTemplate &b) const { return LeftMeet(b); }

    // a.RightMeet(b) returns the right meet of a and b.
    FactorTemplate RightMeet(const FactorTemplate &b) const {
        return FactorTemplate(Underlying.RightMeet(b.Underlying));
    }

    // a.LeftJoin(b) returns the left join of a and b.
    FactorTemplate LeftJoin(const FactorTemplate &b) const {
        return RightComplement()
            .RightMeet(b.RightComplement())
            .LeftComplement();
    }

    // a.RightJoin(b) returns the right join of a and b.
    FactorTemplate RightJoin(const FactorTemplate &b) const {
        return LeftComplement().LeftMeet(b.LeftComplement()).RightComplement();
    }

    // a.IsLeftWeighted(b) returns true if a | b is left weighted, or false
    // otherwise.
    bool IsLeftWeighted(const FactorTemplate &b) const {
        return RightComplement().LeftMeet(b).IsTrivial();
    }

    // a.IsRightWeighted(b) returns true if a | b is right weighted, or false
    // otherwise.
    bool IsRightWeighted(const FactorTemplate &b) const {
        return LeftMeet(b.LeftComplement()).IsTrivial();
    }

    // a.Product(b) returns the product of two factors, under the assumption
    // that it lies below Delta.
    FactorTemplate Product(const FactorTemplate &b) const {
        return FactorTemplate(Underlying.Product(b.Underlying));
    }

    void RightProduct(const FactorTemplate &b) { *this = *this * b; }

    std::size_t Hash() const { return Underlying.Hash(); }

    // a * b is the product of a and b, under the assumption that it lies below
    // Delta. Syntactic sugar for a.Product(b).
    FactorTemplate operator*(const FactorTemplate &b) const { return Product(b); }

    // a.Randomize() sets a to a random factor.
    void Randomize() {

#ifdef RANDOMIZE_ON_ATOMS

        std::vector<FactorTemplate> atoms = Atoms();

        *this = atoms[rand() % int(atoms.size())];

#else

        Underlying.Randomize();

#endif
    }

    // a.Atoms() returns the list of the atoms.
    std::vector<FactorTemplate> Atoms() const {
        std::vector<U> atoms = Underlying.Atoms();
        typename std::vector<U>::iterator atoms_it;
        std::vector<FactorTemplate> factor_atoms;
        for (auto const &atoms_it : atoms) {
            factor_atoms.push_back(atoms_it);
        }
        return factor_atoms;
    }
};

// MakeLeftWeighted(u, v) computes the left-weighted decomposition u' | v' =
// u | v, and sets u = u' and v = v'. It then returns true if something was
// done (so that it may be used with `apply_binfun`). SHOULD NEVER BE CALLED
// UPON u, v IF
// `&u == &v`!
template <class F> bool MakeLeftWeighted(F &u, F &v) {
    F t = (~u) ^ v;
    if (t.IsIdentity()) {
        return false;
    } else {
        v = v / t;
        u = u * t;
        return true;
    }
}

// MakeRightWeighted(u, v) computes the right-weighted decomposition u' | v'
// = u | v, and sets u = u' and v = v'. It then returns true if something
// was done (so that it may be used with `apply_binfun`). SHOULD NEVER BE
// CALLED UPON u, v IF `&u == &v`!
template <class F> bool MakeRightWeighted(F &u, F &v) {
    F t = u.RightMeet(v.LeftComplement());
    if (t.IsIdentity()) {
        return false;
    } else {
        v = t * v;
        u = t.LeftComplement(u);
        return true;
    }
}

// Overloading << for factor classes.
template <class U>
IndentedOStream &operator<<(IndentedOStream &os, const FactorTemplate<U> &f) {
    f.Print(os);
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
    typedef typename F::ParameterType ParameterType;

    /**
     * @brief A (group) parameter.
     *
     * We often implement Garside structures for a series of groups, that
     * are distinguished one from another by some parameter (often, but not
     * always, an integer).
     *
     */
    ParameterType Parameter;

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
    typedef typename std::list<F>::iterator FactorItr;
    typedef typename std::list<F>::const_iterator ConstFactorItr;
    typedef typename std::list<F>::reverse_iterator RevFactorItr;
    typedef typename std::list<F>::const_reverse_iterator ConstRevFactorItr;

  public:
    /**
     * @brief Construct a new BraidTemplate, with a group parameter.
     *
     * Construct a new BraidTemplate, with `Parameter = parameter`. It is
     * initialized as the identity braid.
     *
     * @param parameter Group parameter.
     */
    BraidTemplate(ParameterType parameter)
        : Parameter(parameter), Delta(0), FactorList() {}

    /**
     * @brief Construct a new BraidTemplate, from a factor.
     *
     * Construct a new BraidTemplate whose only factor is `f`.
     *
     * @param f FactorTemplate to be converted to a braid.
     */
    BraidTemplate(const F &f) : Parameter(f.GetParameter()), Delta(0), FactorList() {
        if (f.IsDelta()) {
            Delta = 1;
        } else if (!f.IsIdentity()) {
            FactorList.push_back(f);
        }
    }

    inline static ParameterType parameter_of_string(const std::string &str) {
        return F::parameter_of_string(str);
    }

    ParameterType GetParameter() const { return Parameter; }

    /**
     * @brief Prints `*this` to `os`.
     *
     * Prints `*this` to `os`, in a format that is compatible with
     * `OfString` (assuming that the same holds for `F`).
     *
     * @param os The output stream it prints to.
     */
    void Print(IndentedOStream &os = ind_cout) const {
        if (Delta != 0 && Delta != 1) {
            os << "D ^ " << Delta << (CanonicalLength() > 0 ? " . " : "");
        } else if (Delta == 1) {
            os << "D" << (CanonicalLength() > 0 ? " . " : "");
        }
        ConstFactorItr it_end = FactorList.end();
        it_end--;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            (*it).Print(os);
            if (it != it_end) {
                os << " . ";
            }
        }
    }

    /**
     * @brief Prints `*this` (in RCF) to `os`.
     *
     * Prints `*this` to `os`, in a format that is compatible with
     * `OfString` (assuming that the same holds for `F`).
     *
     * @param os The output stream it prints to.
     */
    void print_rcf(IndentedOStream &os = ind_cout) const {
        ConstFactorItr it_end = FactorList.end();
        it_end--;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            (*it).Print(os);
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

    // `w.Identity()` sets w to the empty word.
    inline void Identity() {
        Delta = 0;
        FactorList.clear();
    }

    // `u.CanonicalLength` returns u's canonical length.
    inline sint16 CanonicalLength() const { return FactorList.size(); }

    inline sint16 Inf() const { return Delta; }

    inline sint16 Sup() const { return Inf() + CanonicalLength(); }

    // `u.Compare(v)` returns whether u and v have the same internal
    // representation.
    inline bool Compare(const BraidTemplate &v) const {
        return (Delta == v.Delta && FactorList == v.FactorList);
    }

    // `u == v` returns whether u and v have the same internal
    // representation. Syntactic sugar for `u.Compare(v)`.
    bool operator==(const BraidTemplate &v) const { return Compare(v); }

    // `u != v` returns whether u and v do not have the same internal
    // representation.
    bool operator!=(const BraidTemplate &v) const { return !Compare(v); }

    // `u.IsIdentity` returns whether u represents the Identity element.
    bool IsIdentity() const { return Delta == 0 && FactorList.empty(); }

    // `u.Inverse()` returns the inverse of u.
    //  See the ElRifai and Morton 1994 article for correction.
    BraidTemplate Inverse() const {
        BraidTemplate b(GetParameter());
        b.Delta = -Delta;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            // Rewrite a_1 ... a_k (f)^(- 1) Delta^r as
            // a_1 ... a_k Delta ^ (r - 1) (Delta^(- r) d_L(f) Delta^r).
            b.FactorList.push_front(
                (*it).LeftComplement().DeltaConjugate(b.Delta));
            --b.Delta;
        }
        return b;
    }

    // `u.Inverse()` returns the inverse of u.
    //  See the ElRifai and Morton 1994 article for correction.
    BraidTemplate InverseRCF() const {
        BraidTemplate b(GetParameter());
        b.Delta = -Delta;
        for (ConstRevFactorItr revit = FactorList.rbegin();
             revit != FactorList.rend(); revit++) {
            // Rewrite Delta^r (f)^(- 1) a_1 ... a_k as
            // (Delta^r d_R(f) Delta^(- r)) Delta ^ (r - 1) a_1 ... a_k.
            b.FactorList.push_back(
                (*revit).RightComplement().DeltaConjugate(-b.Delta));
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
        while (it != FactorList.end() && (*it).IsDelta()) {
            ++it;
            ++Delta;
        }
        FactorList.erase(FactorList.begin(), it);
        RevFactorItr revit = FactorList.rbegin();
        while (revit != FactorList.rend() && (*revit).IsIdentity()) {
            ++revit;
        }
        FactorList.erase(revit.base(), FactorList.end());
    }

    void CleanRCF() {
        FactorItr it = FactorList.begin();
        while (it != FactorList.end() && (*it).IsIdentity()) {
            ++it;
        }
        FactorList.erase(FactorList.begin(), it);
        RevFactorItr revit = FactorList.rbegin();
        while (revit != FactorList.rend() && (*revit).IsDelta()) {
            ++revit;
            ++Delta;
        }
        FactorList.erase(revit.base(), FactorList.end());
    }

    // `u.LeftProduct(f)` assigns fu to u.
    void LeftProduct(const F &f) {
        FactorList.push_front(f.DeltaConjugate(Delta));
        apply_binfun(FactorList.begin(), FactorList.end(), MakeLeftWeighted<F>);
        Clean();
    }

    // `u.RightProduct(f)` assigns uf to u.
    void RightProduct(const F &f) {
        FactorList.push_back(f);
        reverse_apply_binfun(FactorList.begin(), FactorList.end(),
                             MakeLeftWeighted<F>);
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

    // `u.RightProduct(v)` assigns u v to u.
    // v's factors move directly to u - be careful.
    void RightProduct(const BraidTemplate &v) {
        for (FactorItr it = FactorList.begin(); it != FactorList.end(); it++) {
            (*it) = (*it).DeltaConjugate(v.Delta);
        }
        Delta += v.Delta;
        for (ConstFactorItr it = v.FactorList.begin(); it != v.FactorList.end();
             it++) {
            RightProduct((*it));
        }
    }

    // `u.LeftDivide(v)` assigns v ^ (- 1) u to u.
    inline void LeftDivide(const BraidTemplate &v) { LeftProduct(!v); }

    // `u.RightDivide(v)` assigns u v ^ (- 1) to u.
    inline void RightDivide(const BraidTemplate &v) { RightProduct(!v); }

    // `u.LeftDivide(f)` assigns f ^ (- 1) u to u.
    inline void LeftDivide(const F &f) { LeftProduct(!BraidTemplate(f)); }

    // `u.RightDivide(v)` assigns u v ^ (- 1) to u.
    inline void RightDivide(const F &f) { RightProduct(!BraidTemplate(f)); }

    // `u.LeftProductRCF(f)` assigns fu to u.
    void LeftProductRCF(const F &f) {
        FactorList.push_front(f);
        apply_binfun(FactorList.begin(), FactorList.end(),
                     MakeRightWeighted<F>);
        CleanRCF();
    }

    // `u.RightProduct(f)` assigns uf to u.
    void RightProductRCF(const F &f) {
        FactorList.push_back(f.DeltaConjugate(-Delta));
        reverse_apply_binfun(FactorList.begin(), FactorList.end(),
                             MakeRightWeighted<F>);
        CleanRCF();
    }

    // `u.LeftProduct(v)` assigns v u to u.
    void LeftProductRCF(const BraidTemplate &v) {
        for (RevFactorItr it = FactorList.rbegin(); it != FactorList.rend();
             it++) {
            (*it) = (*it).DeltaConjugate(-v.Delta);
        }
        Delta += v.Delta;
        for (ConstRevFactorItr it = v.FactorList.rbegin();
             it != v.FactorList.rend(); it++) {
            LeftProductRCF(*it);
        }
    }

    // `u.RightProduct(v)` assigns u v to u.
    void RightProductRCF(const BraidTemplate &v) {
        for (ConstFactorItr it = v.FactorList.begin(); it != v.FactorList.end();
             it++) {
            RightProductRCF(*it);
        }
        Delta += v.Delta;
    }

    BraidTemplate LeftMeet(const BraidTemplate &v) const {
        sint16 shift = 0;
        BraidTemplate b = BraidTemplate(GetParameter());
        F f1 = F(GetParameter()), f2 = F(GetParameter()), f = F(GetParameter());
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

        while (!f.IsIdentity()) {
            if (b1.Delta > 0) {
                f1.Delta();
            } else if (b1.CanonicalLength() == 0) {
                f1.Identity();
            } else {
                f1 = b1.FactorList.front();
            }

            if (b2.Delta > 0) {
                f2.Delta();
            } else if (b2.CanonicalLength() == 0) {
                f2.Identity();
            } else {
                f2 = b2.FactorList.front();
            }

            f = f1 ^ f2;

            b.RightProduct(f);
            b1.LeftDivide(f);
            b2.LeftDivide(f);
        }

        b.Delta -= shift;
        return b;
    }

    inline BraidTemplate LeftMeet(const F &f) { return LeftMeet(BraidTemplate(f)); }

    inline BraidTemplate operator^(const BraidTemplate &v) { return LeftMeet(v); }

    inline BraidTemplate operator^(const F &f) { return LeftMeet(f); }

    BraidTemplate LeftJoin(const BraidTemplate &v) const {
        sint16 shift = 0;
        BraidTemplate b = BraidTemplate(GetParameter());
        F f2 = F(GetParameter()), f = F(GetParameter());
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

        while (!b2.IsIdentity()) {
            if (b2.Delta > 0) {
                f2.Delta();
            } else if (b2.CanonicalLength() == 0) {
                f2.Identity();
            } else {
                f2 = b2.FactorList.front();
            }

            f = b1.Remainder(f2);

            b.RightProduct(f);
            b1.RightProduct(f);
            b1.LeftDivide(f2);
            b2.LeftDivide(f2);
        }

        b.Delta -= shift;
        return b;
    }

    inline BraidTemplate LeftJoin(const F &f) const { return LeftJoin(BraidTemplate(f)); }

    inline BraidTemplate RightMeet(const BraidTemplate &v) const {
        return !((!(*this)).LeftJoin(!v));
    }

    inline BraidTemplate RightMeet(const F &f) const {
        return !((!(*this)).LeftJoin(!BraidTemplate(f)));
    }

    inline BraidTemplate RightJoin(const BraidTemplate &v) const {
        return !((!(*this)).LeftMeet(!v));
    }

    inline BraidTemplate RightJoin(const F &f) const {
        return !((!(*this)).LeftMeet(!BraidTemplate(f)));
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
        RightProduct(f);
    }

    inline void Conjugate(const BraidTemplate &v) {
        LeftDivide(v);
        RightProduct(v);
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
            F id = F(GetParameter());
            id.Identity();
            return id;
        } else {
            return FactorList.front().DeltaConjugate(-Delta);
        }
    }

    inline F First() const {
        if (CanonicalLength() == 0) {
            F id = F(GetParameter());
            id.Identity();
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
            F id = F(GetParameter());
            id.Identity();
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
            F id = F(GetParameter());
            id.Identity();
            return id;
        } else {
            return FactorList.back().DeltaConjugate(Delta).RightMeet(
                FactorList.front().LeftComplement());
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
        RightProduct(i);
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

    // `u.Product(v)` returns uv.
    BraidTemplate Product(const BraidTemplate &v) const {
        BraidTemplate w(*this);
        w.RightProduct(v);
        return w;
    }

    // `u.Power(k)` returns u raised to the power k.
    // Uses a fast exponentiation algorithm; the number of multiplications
    // is logarithmic in k.
    BraidTemplate Power(const sint16 k) const {
        if (k == 0) {
            return BraidTemplate(GetParameter());
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
    // Syntactic sugar for `u.Product(v)`.
    BraidTemplate operator*(const BraidTemplate &v) const { return Product(v); }

    // `u.Normalize()` turns u into LCF.
    inline void Normalize() {
        bubble_sort(FactorList.begin(), FactorList.end(), MakeLeftWeighted<F>);
        Clean();
    }

    // `u.LCFToRCF()` turns u, assumed to be in LCF (so that we may avoid
    // having to clean up), into RCF. The result (r, [u_1, ..., u_k])
    // represents u_1
    // ... u_k Delta ^ r.
    inline void MakeRCFFromLCF() {
        for (FactorItr it = FactorList.begin(); it != FactorList.end(); ++it) {
            *it = (*it).DeltaConjugate(-Delta);
        };
        bubble_sort(FactorList.begin(), FactorList.end(), MakeRightWeighted<F>);
    }

    // `u.LCFToRCF()` turns u, assumed to be in LCF, into RCF (so that we
    // may avoid having to clean up). The result (r, [u_1, ..., u_k])
    // represents u_1
    // ... u_k Delta ^ r.
    inline void MakeLCFFromRCF() {
        for (FactorItr it = FactorList.begin(); it != FactorList.end(); ++it) {
            *it = (*it).DeltaConjugate(Delta);
        };
        bubble_sort(FactorList.begin(), FactorList.end(), MakeLeftWeighted<F>);
    }

    // `b.Remainder(f)` computes, if b is positive, the simple factor s such
    // that bs is the left lcm of b and b and f.
    F Remainder(const F &f) const {
        F fi = f;
        if (Delta != 0) {
            fi.Identity();
        } else {
            for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
                 it++) {
                fi = (*it).LeftJoin(fi) / *it;
            }
        }
        return fi;
    }

    sint16 Rigidity() const {
        BraidTemplate b2 = *this;
        sint16 rigidity = 0;

        if (CanonicalLength() == 0)
            return rigidity;

        b2.RightProduct(b2.Initial());

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
    void Randomize(sint16 canonical_length) {
        FactorList.clear(); // Potential memory leak? To be checked.
        Delta = 0;
        for (sint16 i = 0; i < canonical_length; i++) {
            F f = F(GetParameter());
            f.Randomize();
            FactorList.push_back(f);
        }
    }

    void Debug(IndentedOStream &os = ind_cout) const {
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
            (*it).Debug(os);
            if (it != it_end) {
                os << "," << EndLine();
            }
        }
        os.Indent(-4);
        os << EndLine() << "]";
        os.Indent(-8);
        os << EndLine() << "}";
    }

    std::size_t Hash() const {
        std::size_t h = Delta;
        for (ConstFactorItr it = FactorList.begin(); it != FactorList.end();
             it++) {
            h = h * 31 + (*it).Hash();
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
    void OfString(const std::string str) {
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
        b.Identity();

        std::smatch match;

        F fact(GetParameter());

        std::regex_search(str.begin() + pos, str.end(), match, ignore,
                          std::regex_constants::match_continuous);
        pos += match[0].length();

        while (pos != str.length()) {
            sint16 pow;

            fact.OfString(str, pos);

            if (std::regex_search(str.begin() + pos, str.end(), match, power,
                                  std::regex_constants::match_continuous)) {
                pos += match[0].length();
                pow = std::stoi(match[1]);
            } else {
                pow = 1;
            }
            if (pow >= 0) {
                for (sint16 _ = 0; _ < pow; _++) {
                    b.RightProduct(fact);
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
    b.Print(os);
    return os;
}

} // namespace cgarside

template <class U> struct std::hash<cgarside::FactorTemplate<U>> {
    std::size_t operator()(cgarside::FactorTemplate<U> const &f) const noexcept {
        return f.Hash();
    }
};

template <class F> struct std::hash<cgarside::BraidTemplate<F>> {
    std::size_t operator()(cgarside::BraidTemplate<F> const &u) const noexcept {
        return u.Hash();
    }
};

#endif
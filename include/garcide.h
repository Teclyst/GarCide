/**
 * @file garcide.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header (and implementation) file for generic Garside groups.
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

#ifndef GARCIDE
#define GARCIDE

#include "utility.h"
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace garcide {

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
        return left_complement()
            .left_meet(b.left_complement())
            .right_complement();
    }

    // a.is_left_weighted(b) returns true if a | b is left weighted, or false
    // otherwise.
    bool is_left_weighted(const FactorTemplate &b) const {
        return right_complement().left_meet(b).is_identity();
    }

    // a.is_right_weighted(b) returns true if a | b is right weighted, or false
    // otherwise.
    bool is_right_weighted(const FactorTemplate &b) const {
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

  private:
    /**
     * @brief A (group) parameter.
     *
     * We often implement Garside structures for a series of groups, that
     * are distinguished one from another by some parameter (often, but not
     * always, an integer).
     */
    Parameter parameter;
    /**
     * @brief Infimum.
     *
     * The power of delta at the left end of the word (in LCF, otherwise at
     * the right end in RCF).
     *
     */
    sint32 delta;

    /**
     * @brief The braid's canonical factors.
     *
     * A list of the braid's canonical factors, from left to right. It is
     * left weighted when in LCF, and right weighted when in RCF.
     */
    std::list<F> factor_list;

  public:
    using FactorItr = typename std::list<F>::iterator;
    using RevFactorItr = typename std::list<F>::reverse_iterator;
    using ConstFactorItr = typename std::list<F>::const_iterator;
    using ConstRevFactorItr = typename std::list<F>::const_reverse_iterator;

    inline FactorItr begin() { return factor_list.begin(); }

    inline RevFactorItr rbegin() { return factor_list.rbegin(); }

    inline ConstFactorItr cbegin() const { return factor_list.begin(); }

    inline ConstRevFactorItr crbegin() const { return factor_list.rbegin(); }

    inline FactorItr end() { return factor_list.end(); }

    inline RevFactorItr rend() { return factor_list.rend(); }

    inline ConstFactorItr cend() const { return factor_list.end(); }

    inline ConstRevFactorItr crend() const { return factor_list.rend(); }

  public:
    /**
     * @brief Construct a new BraidTemplate, with a group parameter.
     *
     * Construct a new BraidTemplate, with `parameter = parameter`. It is
     * initialized as the identity braid.
     *
     * @param parameter Group parameter.
     */
    BraidTemplate(Parameter parameter)
        : parameter(parameter), delta(0), factor_list() {}

    /**
     * @brief Construct a new BraidTemplate, from a factor.
     *
     * Construct a new BraidTemplate whose only factor is `f`.
     *
     * @param f FactorTemplate to be converted to a braid.
     */
    BraidTemplate(const F &f)
        : parameter(f.get_parameter()), delta(0), factor_list() {
        if (f.is_delta()) {
            delta = 1;
        } else if (!f.is_identity()) {
            factor_list.push_back(f);
        }
    }

    inline static Parameter parameter_of_string(const std::string &str) {
        return F::parameter_of_string(str);
    }

    Parameter get_parameter() const { return parameter; }

    /**
     * @brief Sets `delta` to `delta`.
     *
     * Sets `delta` field (which is private) to `delta`.
     *
     * @param delta The new value of `delta`.
     */
    inline void set_delta(sint16 delta) { (*this).delta = delta; }

    /**
     * @brief Prints `*this` to `os`.
     *
     * Prints `*this` to `os`, in a format that is compatible with
     * `of_string` (assuming that the same holds for `F`).
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os = ind_cout) const {
        if (delta != 0 && delta != 1) {
            os << "D ^ " << delta << (canonical_length() > 0 ? " . " : "");
        } else if (delta == 1) {
            os << "D" << (canonical_length() > 0 ? " . " : "");
        }
        ConstFactorItr it_end = cend();
        it_end--;
        for (ConstFactorItr it = cbegin(); it != cend(); it++) {
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
        ConstFactorItr it_end = cend();
        it_end--;
        for (ConstFactorItr it = cbegin(); it != cend(); it++) {
            (*it).print(os);
            if (it != it_end) {
                os << " . ";
            }
        }
        if (delta != 0 && delta != 1) {
            os << (canonical_length() > 0 ? " . " : "") << "D ^ " << delta;
        } else if (delta == 1) {
            os << (canonical_length() > 0 ? " . " : "") << "D";
        }
    }

    // `w.identity()` sets w to the empty word.
    inline void identity() {
        delta = 0;
        factor_list.clear();
    }

    // `u.canonical_length` returns u's canonical length.
    inline size_t canonical_length() const { return factor_list.size(); }

    inline sint32 inf() const { return delta; }

    inline sint32 sup() const { return inf() + int(canonical_length()); }

    // `u.compare(v)` returns whether u and v have the same internal
    // representation.
    inline bool compare(const BraidTemplate &v) const {
        return (delta == v.delta && factor_list == v.factor_list);
    }

    // `u == v` returns whether u and v have the same internal
    // representation. Syntactic sugar for `u.compare(v)`.
    bool operator==(const BraidTemplate &v) const { return compare(v); }

    // `u != v` returns whether u and v do not have the same internal
    // representation.
    bool operator!=(const BraidTemplate &v) const { return !compare(v); }

    // `u.is_identity` returns whether u represents the identity element.
    bool is_identity() const { return delta == 0 && factor_list.empty(); }

    // `u.inverse()` returns the inverse of u.
    //  See the ElRifai and Morton 1994 article for correction.
    BraidTemplate inverse() const {
        BraidTemplate b(get_parameter());
        b.delta = -delta;
        for (ConstFactorItr it = cbegin(); it != cend(); it++) {
            // Rewrite a_1 ... a_k (f)^(- 1) Delta^r as
            // a_1 ... a_k Delta ^ (r - 1) (Delta^(- r) d_L(f) Delta^r).
            b.factor_list.push_front(
                (*it).left_complement().delta_conjugate(b.delta));
            --b.delta;
        }
        return b;
    }

    // `u.inverse()` returns the inverse of u.
    //  See the ElRifai and Morton 1994 article for correction.
    BraidTemplate inverse_rcf() const {
        BraidTemplate b(get_parameter());
        b.delta = -delta;
        for (ConstRevFactorItr revit = crbegin(); revit != crend(); revit++) {
            // Rewrite Delta^r (f)^(- 1) a_1 ... a_k as
            // (Delta^r d_R(f) Delta^(- r)) Delta ^ (r - 1) a_1 ... a_k.
            b.factor_list.push_back(
                (*revit).right_complement().delta_conjugate(-b.delta));
            --b.delta;
        }
        return b;
    }

    // `!u` returns the inverse of u.
    // Syntactic sugar for `u.inverse()`.
    inline BraidTemplate operator!() const { return inverse(); }

    // clean gets rid of (factor) Deltas at the start, and identity elements
    // at the end pf `factor_list`.
    void clean() {
        FactorItr it = begin();
        while (it != end() && (*it).is_delta()) {
            ++it;
            ++delta;
        }
        factor_list.erase(begin(), it);
        RevFactorItr revit = rbegin();
        while (revit != rend() && (*revit).is_identity()) {
            ++revit;
        }
        factor_list.erase(revit.base(), end());
    }

    void clean_rcf() {
        FactorItr it = begin();
        while (it != end() && (*it).is_identity()) {
            ++it;
        }
        factor_list.erase(begin(), it);
        RevFactorItr revit = rbegin();
        while (revit != rend() && (*revit).is_delta()) {
            ++revit;
            ++delta;
        }
        factor_list.erase(revit.base(), end());
    }

    // `u.left_multiply(f)` assigns fu to u.
    inline void left_multiply(const F &f) {
        factor_list.push_front(f.delta_conjugate(delta));
        apply_binfun(begin(), end(), make_left_weighted<F>);
        clean();
    }

    // `u.right_multiply(f)` assigns uf to u.
    inline void right_multiply(const F &f) {
        factor_list.push_back(f);
        reverse_apply_binfun(begin(), end(), make_left_weighted<F>);
        clean();
    }

    // `u.left_multiply(v)` assigns v u to u.
    void left_multiply(const BraidTemplate &v) {
        for (ConstRevFactorItr it = v.crbegin(); it != v.crend(); it++) {
            left_multiply(*it);
        }
        delta += v.delta;
    }

    // `u.right_multiply(v)` assigns u v to u.
    void right_multiply(const BraidTemplate &v) {
        for (FactorItr it = begin(); it != end(); it++) {
            (*it).delta_conjugate_mut(v.delta);
        }
        delta += v.delta;
        for (ConstFactorItr it = v.cbegin(); it != v.cend(); it++) {
            right_multiply((*it));
        }
    }

    // `u.left_divide(v)` assigns v ^ (- 1) u to u.
    inline void left_divide(const BraidTemplate &v) { left_multiply(!v); }

    // `u.right_divide(v)` assigns u v ^ (- 1) to u.
    inline void right_divide(const BraidTemplate &v) { right_multiply(!v); }

    // `u.left_divide(f)` assigns f ^ (- 1) u to u.
    inline void left_divide(const F &f) { left_multiply(!BraidTemplate(f)); }

    // `u.right_divide(v)` assigns u v ^ (- 1) to u.
    inline void right_divide(const F &f) { right_multiply(!BraidTemplate(f)); }

    // `u.left_multiply_rcf(f)` assigns fu to u.
    inline void left_multiply_rcf(const F &f) {
        factor_list.push_front(f);
        apply_binfun(begin(), end(), make_right_weighted<F>);
        clean_rcf();
    }

    // `u.right_multiply(f)` assigns uf to u.
    inline void right_multiply_rcf(const F &f) {
        factor_list.push_back(f.delta_conjugate(-delta));
        reverse_apply_binfun(begin(), end(), make_right_weighted<F>);
        clean_rcf();
    }

    // `u.left_multiply(v)` assigns v u to u.
    void left_multiply_rcf(const BraidTemplate &v) {
        for (RevFactorItr it = rbegin(); it != rend(); it++) {
            (*it).delta_conjugate_mut(-v.delta);
        }
        delta += v.delta;
        for (ConstRevFactorItr it = v.crbegin(); it != v.crend(); it++) {
            left_multiply_rcf(*it);
        }
    }

    // `u.right_multiply(v)` assigns u v to u.
    void right_multiply_rcf(const BraidTemplate &v) {
        for (ConstFactorItr it = v.cbegin(); it != v.cend(); it++) {
            right_multiply_rcf(*it);
        }
        delta += v.delta;
    }

    // `u.left_divide(v)` assigns v ^ (- 1) u to u.
    inline void left_divide_rcf(const BraidTemplate &v) {
        left_multiply_rcf(v.inverse_rcf());
    }

    // `u.right_divide(v)` assigns u v ^ (- 1) to u.
    inline void right_divide_rcf(const BraidTemplate &v) {
        right_multiply_rcf(v.inverse_rcf());
    }

    // `u.left_divide(f)` assigns f ^ (- 1) u to u.
    inline void left_divide_rcf(const F &f) {
        left_multiply_rcf(BraidTemplate(f).inverse_rcf());
    }

    // `u.right_divide(v)` assigns u v ^ (- 1) to u.
    inline void right_divide_rcf(const F &f) {
        right_multiply_rcf(BraidTemplate(f).inverse_rcf());
    }

    BraidTemplate left_meet(const BraidTemplate &v) const {
        sint16 shift = 0;
        BraidTemplate b(get_parameter());
        F f1(get_parameter()), f2(get_parameter()), f(get_parameter());
        f.delta();

        BraidTemplate b1 = *this, b2 = v;

        shift -= b1.delta;
        b2.delta -= b1.delta;
        b1.delta = 0;

        if (b2.delta < 0) {
            shift -= b2.delta;
            b1.delta -= b2.delta;
            b2.delta = 0;
        }

        while (!f.is_identity()) {
            if (b1.delta > 0) {
                f1.delta();
            } else if (b1.canonical_length() == 0) {
                f1.identity();
            } else {
                f1 = b1.first();
            }

            if (b2.delta > 0) {
                f2.delta();
            } else if (b2.canonical_length() == 0) {
                f2.identity();
            } else {
                f2 = b2.first();
            }

            f = f1 ^ f2;

            b.right_multiply(f);
            b1.left_divide(f);
            b2.left_divide(f);
        }

        b.delta -= shift;
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
        f.delta();

        BraidTemplate b1 = *this, b2 = v;

        shift -= b1.delta;
        b2.delta -= b1.delta;
        b1.delta = 0;

        if (b2.delta < 0) {
            shift -= b2.delta;
            b1.delta -= b2.delta;
            b2.delta = 0;
        }

        b = b1;

        while (!b2.is_identity()) {
            if (b2.inf() > 0) {
                f2.delta();
            } else if (b2.canonical_length() == 0) {
                f2.identity();
            } else {
                f2 = b2.first();
            }

            f = b1.remainder(f2);

            b.right_multiply(f);
            b1.right_multiply(f);
            b1.left_divide(f2);
            b2.left_divide(f2);
        }

        b.delta -= shift;
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

    inline void conjugate(const F &f) {
        left_divide(f);
        right_multiply(f);
    }

    inline void conjugate(const BraidTemplate &v) {
        left_divide(v);
        right_multiply(v);
    }

    inline void conjugate_rcf(const F &f) {
        left_divide_rcf(f);
        right_multiply_rcf(f);
    }

    inline void conjugate_rcf(const BraidTemplate &v) {
        left_divide_rcf(v);
        right_multiply_rcf(v);
    }

    /**
     * @brief Returns the first (non-Delta) factor.
     *
     * Returns the first (non-Delta) factor. If canonical length is zero,
     * returns the identity factor instead.
     *
     * @return The first (non-Delta) factor.
     */
    inline F first() const {
        if (canonical_length() == 0) {
            F id(get_parameter());
            id.identity();
            return id;
        } else {
            return factor_list.front();
        }
    }

    /**
     * @brief Returns the initial factor.
     *
     * Returns the initial factor, which is the first factor conjugated by
     * `-inf()`. If canonical length is zero, returns the identity factor
     * instead.
     *
     * @return The first (non-Delta) factor.
     */
    inline F initial() const { return first().delta_conjugate(-inf()); }

    /**
     * @brief Returns the final factor.
     *
     * Returns the final factor, (i.e. the last one in `factor_list`). If
     * canonical length is zero, returns the identity factor instead.
     *
     * @return The first (non-Delta) factor.
     */
    inline F final() const {
        if (canonical_length() == 0) {
            F id(get_parameter());
            id.identity();
            return id;
        } else {
            return factor_list.back();
        }
    }

    // `u.preferred_prefix()` returns the preferred prefix of u, that is, if
    // u = Delta ^ r u_1 ... u_k, p(u) = d_R(u_k) ^_L Delta ^ r u_1 Delta ^
    // (- r) If u has canonical length zero, returns the identity factor
    // instead.
    inline F preferred_prefix() const { return initial() ^ ~final(); }

    F preferred_suffix_rcf() const {
        if (canonical_length() == 0) {
            F id(get_parameter());
            id.identity();
            return id;
        } else {
            return factor_list.back().delta_conjugate(inf()).right_meet(
                factor_list.front().left_complement());
        }
    };

    inline F preferred_suffix() const {
        BraidTemplate right = *this;
        right.lcf_to_rcf();
        return right.preferred_suffix_rcf();
    };

    // `u.cycling()` cycles u: if u = Delta ^ r u_1 ... u_k, then after
    // applying cycling u will contain (the LNF of) Delta ^ r u_2 ... u_k
    // (Delta ^ r u_1 Delta ^ (-r)).
    inline void cycling() {
        if (canonical_length() == 0) {
            return;
        }
        F i = initial();
        factor_list.pop_front();
        right_multiply(i);
    }

    // `u.decycling()` decycles u: if u = Delta ^ r u_1 ... u_k, then after
    // applying cycling u will contain (the LNF of) Delta ^ r u_2 ... u_k
    // (Delta ^ r u_1 Delta ^ (-r)).
    inline void decycling() {
        if (canonical_length() == 0) {
            return;
        }
        F f = final();
        factor_list.pop_back();
        left_multiply(f);
    }

    // `u.sliding()` cyclically slides u: if u = Delta ^ r u_1 ... u_k, and
    // Delta ^ r u_1 Delta ^ (-r) = p(u) u'_1, then after applying cycling u
    // will contain (the LNF of) Delta ^ r (Delta ^ r u'_1 Delta ^ (-r)) u_2
    // ... u_k p(u).
    inline void sliding() {
        if (canonical_length() == 0) {
            return;
        }
        conjugate(preferred_prefix());
    }

    // `u.product(v)` returns uv.
    BraidTemplate product(const BraidTemplate &v) const {
        BraidTemplate w(*this);
        w.right_multiply(v);
        return w;
    }

    // `u.power(k)` returns u raised to the power k.
    // Uses a fast exponentiation algorithm; the number of multiplications
    // is logarithmic in k.
    /**
     * @brief Raises to `k`-th power.
     *
     * Computes `*this` raised to the `k`-th power, using a fast exponentiation
     * algorithm.
     *
     * @param k The exponent.
     * @return `*this` raised to the `k`-th power.
     */
    BraidTemplate power(const sint16 k) const {
        if (k == 0) {
            return BraidTemplate(get_parameter());
        } else if (k % 2 == 0) {
            BraidTemplate root = power(k / 2);
            return root * root;
        } else if (k > 0) {
            BraidTemplate root = power(k / 2);
            return *this * root * root;
        } else {
            BraidTemplate root = power(k / 2);
            return !*this * root * root;
        }
    }

    // `u * v` returns uv.
    // Syntactic sugar for `u.product(v)`.
    BraidTemplate operator*(const BraidTemplate &v) const { return product(v); }

    // `u.normalize()` turns u into LCF.
    inline void normalize() {
        bubble_sort(begin(), end(), make_left_weighted<F>);
        clean();
    }

    /**
     * @brief Changes orientation to RCF.
     *
     * Changes convention: Deltas are now on the right. Also (right-)normalizes
     * the braid.
     */
    inline void lcf_to_rcf() {
        for (FactorItr it = begin(); it != end(); ++it) {
            *it = (*it).delta_conjugate(-delta);
        };
        bubble_sort(begin(), end(), make_right_weighted<F>);
    }

    /**
     * @brief Changes orientation to LCF.
     *
     * Changes convention: Deltas are now on the left. Also (left-)normalizes
     * the braid.
     */
    inline void rcf_to_lcf() {
        for (FactorItr it = begin(); it != end(); ++it) {
            *it = (*it).delta_conjugate(delta);
        };
        bubble_sort(begin(), end(), make_left_weighted<F>);
    }

    // `b.remainder(f)` computes, if b is positive, the simple factor s such
    // that bs is the left lcm of b and b and f.
    F remainder(const F &f) const {
        F fi = f;
        if (delta != 0) {
            fi.identity();
        } else {
            for (ConstFactorItr it = cbegin(); it != cend(); it++) {
                fi = (*it).left_join(fi) / *it;
            }
        }
        return fi;
    }

    sint16 rigidity() const {
        BraidTemplate b2 = *this;
        sint16 rigidity = 0;

        if (canonical_length() == 0)
            return rigidity;

        b2.right_multiply(b2.initial());

        ConstFactorItr it2 = b2.cbegin();

        for (ConstFactorItr it = cbegin(); it != cend();
             it++, it2++, rigidity++) {
            if (*it != *it2) {
                break;
            }
        }

        return rigidity;
    }

    // Randomizes the braid, setting it at a given length. The result isn't in
    // LCF (so that it may be used to benchmark `normalize`).
    /**
     * @brief Randomizes the braid.
     *
     * Randomizes the braid, setting it at a given length. The result isn't in
     * LCF (so that it may be used to benchmark `normalize`).
     * @param canonical_length An upper bound on the braid's new canonical
     * length (number of pushed factors).
     */
    void randomize(sint16 canonical_length) {
        identity();
        F f(get_parameter());
        for (sint16 i = 0; i < canonical_length; i++) {
            f.randomize();
            factor_list.push_back(f);
        }
    }

    void debug(IndentedOStream &os = ind_cout) const {
        os << "{   ";
        os.Indent(4);
        os << "parameter:";
        os.Indent(4);
        os << EndLine();
        os << parameter;
        os.Indent(-4);
        os << EndLine();
        os << "delta:";
        os.Indent(4);
        os << EndLine();
        os << delta;
        os.Indent(-4);
        os << EndLine();
        os << "factor_list:";
        os.Indent(4);
        os << EndLine();
        os << "[   ";
        os.Indent(4);
        ConstFactorItr it_end = cend();
        it_end--;
        for (ConstFactorItr it = cbegin(); it != cend(); it++) {
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
        std::size_t h = inf();
        for (ConstFactorItr it = cbegin(); it != cend(); it++) {
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
                    b.right_divide(fact);
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
inline IndentedOStream &operator<<(IndentedOStream &os,
                                   const BraidTemplate<F> &b) {
    b.print(os);
    return os;
}

} // namespace garcide

template <class U> struct std::hash<garcide::FactorTemplate<U>> {
    std::size_t
    operator()(garcide::FactorTemplate<U> const &f) const noexcept {
        return f.hash();
    }
};

template <class F> struct std::hash<garcide::BraidTemplate<F>> {
    std::size_t operator()(garcide::BraidTemplate<F> const &u) const noexcept {
        return u.hash();
    }
};

#endif
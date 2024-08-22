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

#include "garcide/utility.hpp"
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Namespace for the _GarCide_ library.
 */
namespace garcide {

/**
 * @brief A class template for Garside group canonical factors.
 *
 * @tparam U A template class representing the internal data structure of
 * factors.
 */
template <class U> class FactorTemplate {

  public:
    /**
     * @brief Parameter type.
     *
     * Used to discriminate between the different
     * groups in the family that is being implemented.
     */
    using Parameter = typename U::Parameter;

  private:
    /**
     * @brief The actual data structure representing the factor.
     *
     * `*this` is basically a wrapper for `underlying` (but with more member
     * functions).
     */
    U underlying;

  public:
    /**
     * @brief Construct a new `FactorTemplate`, from an underlying object.
     *
     * The resulting factor is a wrapper for `under`.
     *
     * @param under The underlying object this factor should wrap around.
     */
    FactorTemplate(const U &under) : underlying(under) {}

    /**
     * @brief Construct a new `FactorTemplate`, from a `Parameter` value.
     *
     * The `Parameter` `parameter` will be passed to the corresponding `U`
     * constructor.
     *
     * @param parameter The parameter the `factorTemplate` should have.
     */
    FactorTemplate(Parameter parameter) : underlying(parameter) {}

    /**
     * @brief Gets the underlying representation of the factor.
     *
     * @return `this.underlying`.
     */
    inline U get_underlying() const { return underlying; }

    /**
     * @brief Converts a string to a `Parameter`.
     *
     * This is a wrapper for the matching `U` static member function.
     *
     * @param str The string to be converted.
     * @return A parameter matching `str`.
     * @exception InvalidStringError Raised if `str` cannot be converted to a
     * `Parameter`.
     */
    inline static Parameter parameter_of_string(const std::string &str) {
        return U::parameter_of_string(str);
    }

    /**
     * @brief Returns the `Parameter` of the factor.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @return The `Parameter` of the factor.
     */
    inline Parameter get_parameter() const {
        return underlying.get_parameter();
    }

    /**
     * @brief Computes the height of the factors lattice.
     *
     * That is to say, the maximum length of the Garside element as a product of
     * atoms.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @return The height of the factors lattice.
     */
    inline sint16 lattice_height() const {
        return underlying.lattice_height();
    };

    /**
     * @brief Sets this factor to one matching a substring.
     *
     * Tries to extract a substring matching a factor from `str`, starting at
     * `pos`. If successful, `this` is set to the matching factor and `pos` is
     * advanced by the length of the substring.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param str The string a factor is to be extracted from.
     * @param pos The position from which a substring should be searched.
     *
     * @exception InvalidStringError Thrown when no factor can be extracted
     * starting at `pos`, or a factor can be extracted, but is not legal.
     */
    inline void of_string(const std::string &str, size_t &pos) {
        underlying.of_string(str, pos);
    }

    /**
     * @brief Prints the internal structure of the factor to output stream `os`.
     *
     * Typically used while debugging.
     *
     * @param os The `IndentedOStream` the factor should be printed in.
     */
    void debug(IndentedOStream &os = ind_cout) const {
        os << "{   Underlying:";
        os.Indent(8);
        os << EndLine();
        underlying.debug(os);
        os.Indent(-8);
        os << EndLine() << "}";
    }

    /**
     * @brief Prints the factor to output stream `os`.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param os The `IndentedOStream` the factor should be printed in.
     */
    inline void print(IndentedOStream &os = ind_cout) const {
        underlying.print(os);
    }

    /**
     * @brief Sets the factor to the identity.
     *
     * This is a wrapper for the matching `U` member function.
     */
    inline void identity() { underlying.identity(); }

    /**
     * @brief Sets the factor to the Garside element.
     *
     * This is a wrapper for the matching `U` member function.
     */
    inline void delta() { underlying.delta(); }

    /**
     * @brief Equality test between two factors.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param b The factor `*this` is compared to.
     * @return if `*this` and `b` are equal.
     */
    inline bool compare(const FactorTemplate &b) const {
        return underlying.compare(b.underlying);
    }

    /**
     * @brief Equality test between two factors.
     *
     * Syntactic sugar for `compare`.
     *
     * @param b The factor `*this` is compared to.
     * @return if `*this` and `b` are equal.
     */
    inline bool operator==(const FactorTemplate &b) const { return compare(b); }

    /**
     * @brief Unequality test between two factors.
     *
     * @param b The factor `*this` is compared to.
     * @return if `*this` and `b` are not equal.
     */
    inline bool operator!=(const FactorTemplate &b) const {
        return !compare(b);
    }

    /**
     * @brief Equality test with the identity.
     *
     * @return if `*this` is the identity.
     */
    inline bool is_identity() const {
        FactorTemplate e = FactorTemplate(*this);
        e.identity();
        return compare(e);
    }

    /**
     * @brief Equality test with the Garside element.
     *
     * @return if `*this` is the Garside element.
     */
    inline bool is_delta() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.delta();
        return compare(delta);
    }

    /**
     * @brief Computes the left complement of `*this` under `b`.
     *
     * It is assumed that `*this` right-divides `b`.
     *
     * Given two factors \f$a\f$ and \f$b\f$, with \f$a\f$ a right-divisor of
     * \f$b\f$, the left complement of \f$a\f$ under \f$b\f$ is the factor
     * \f$c\f$ such that \f$b=ca\f$.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param b The factor under which we want to compute the complement.
     * @return The left complement of `*this` under `b`.
     */
    inline FactorTemplate left_complement(const FactorTemplate &b) const {
        return FactorTemplate(underlying.left_complement(b.underlying));
    }

    /**
     * @brief Computes the left complement of `*this` under the Garside element.
     *
     * @return The left complement of `*this` under the Garside element.
     */
    inline FactorTemplate left_complement() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.delta();
        return left_complement(delta);
    }

    /**
     * @brief Computes the right complement of `*this` under `b`.
     *
     * It is assumed that `*this` left-divides `b`.
     *
     * Given two factors \f$a\f$ and \f$b\f$, with \f$a\f$ a left-divisor of
     * \f$b\f$, the right complement of \f$a\f$ under \f$b\f$ is the factor
     * \f$c\f$ such that \f$b=ac\f$.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param b The factor under which we want to compute the complement.
     * @return The right complement of `*this` under `b`.
     */
    inline FactorTemplate right_complement(const FactorTemplate &b) const {
        return FactorTemplate(underlying.right_complement(b.underlying));
    }

    /**
     * @brief Computes the right complement of `*this` under the Garside
     * element.
     *
     * @return The right complement of `*this` under the Garside element.
     */
    inline FactorTemplate right_complement() const {
        FactorTemplate delta = FactorTemplate(*this);
        delta.delta();
        return right_complement(delta);
    }

    /**
     * @brief Computes the right complement of `*this` under the Garside
     * element.
     *
     * Syntactic sugar for `right_complement` (the unary version).
     *
     * @return The right complement of `*this` under the Garside element.
     */
    inline FactorTemplate operator~() const { return right_complement(); }

    /**
     * @brief Computes the right complement of `*this` under `b`.
     *
     * It is assumed that `*this` left-divides `b`.
     *
     * Syntactic sugar for `right_complement` (the binary version).
     *
     * @param b The factor under which we want to compute the complement.
     * @return The right complement of `*this` under `b`.
     */
    inline FactorTemplate operator/(const FactorTemplate &b) const {
        return b.right_complement(*this);
    }

    /**
     * @brief Conjugates the factor by a power of the Garside element, modifying
     * it.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param k The power of the Garside element that the factor should be
     * conjugated by.
     */
    inline void delta_conjugate_mut(sint16 k) {
        underlying.delta_conjugate_mut(k);
    }

    /**
     * @brief Computes the conjugate of the factor by a power of the Garside
     * element.
     *
     * By default, conjugates by the Garside element.
     *
     * @param k The power of the Garside element that the factor should be
     * conjugated by.
     * @return The conjugate.
     */
    FactorTemplate delta_conjugate(sint16 k = 1) const {
        FactorTemplate conjugate = *this;
        conjugate.delta_conjugate_mut(k);
        return conjugate;
    }

    /**
     * @brief Computes left meets.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param b Second operand.
     * @return The left meet of `*this` and `b`.
     */
    inline FactorTemplate left_meet(const FactorTemplate &b) const {
        return FactorTemplate(underlying.left_meet(b.underlying));
    }

    /**
     * @brief Computes left meets.
     *
     * Syntactic sugar for `left_meet`.
     *
     * @param b Second operand.
     * @return The left meet of `*this` and `b`.
     */
    inline FactorTemplate operator^(const FactorTemplate &b) const {
        return left_meet(b);
    }

    /**
     * @brief Computes right meets.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @param b Second operand.
     * @return The right meet of `*this` and `b`.
     */
    inline FactorTemplate right_meet(const FactorTemplate &b) const {
        return FactorTemplate(underlying.right_meet(b.underlying));
    }

    /**
     * @brief Computes left joins.
     *
     * @param b Second operand.
     * @return The left join of `*this` and `b`.
     */
    inline FactorTemplate left_join(const FactorTemplate &b) const {
        return right_complement()
            .right_meet(b.right_complement())
            .left_complement();
    }

    /**
     * @brief Computes right joins.
     *
     * @param b Second operand.
     * @return The right join of `*this` and `b`.
     */
    inline FactorTemplate right_join(const FactorTemplate &b) const {
        return left_complement()
            .left_meet(b.left_complement())
            .right_complement();
    }

    /**
     * @brief Checks left-weightedness.
     *
     * @param b The second factor in the decomposition.
     * @return If `*this`\f${}\mid{}\f$`b` is left-weighted.
     */
    inline bool is_left_weighted(const FactorTemplate &b) const {
        return right_complement().left_meet(b).is_identity();
    }

    /**
     * @brief Checks right-weightedness.
     *
     * @param b The second factor in the decomposition.
     * @return If `*this`\f${}\mid{}\f$`b` is right-weighted.
     */
    inline bool is_right_weighted(const FactorTemplate &b) const {
        return left_meet(b.left_complement()).is_identity();
    }

    /**
     * @brief Computes products.
     *
     * It is assumed that the product is still a factor.
     *
     * @param b Second (right) operand.
     * @return The product of `*this` and `b`.
     */
    inline FactorTemplate product(const FactorTemplate &b) const {
        return FactorTemplate(underlying.product(b.underlying));
    }

    /**
     * @brief Mutably left-multiply.
     *
     * It is assumed that the product is still a factor.
     *
     * @param b The factor we left-multiply by.
     */
    inline void left_multiply(const FactorTemplate &b) { *this = b * *this; }

    /**
     * @brief Mutably right-multiply.
     *
     * It is assumed that the product is still a factor.
     *
     * @param b The factor we right-multiply by.
     */
    inline void right_multiply(const FactorTemplate &b) { *this = *this * b; }

    /**
     * @brief Computes products.
     *
     * It is assumed that the product is still a factor.
     *
     * Syntactic sugar for `product`.
     *
     * @param b Second (right) operand
     * @return The product of `*this` and `b`.
     */
    inline FactorTemplate operator*(const FactorTemplate &b) const {
        return product(b);
    }

    /**
     * @brief Hashes the factor.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @return The hash.
     */
    inline std::size_t hash() const { return underlying.hash(); }

    /**
     * @brief Randomizes the factor.
     *
     * This is a wrapper for the matching `U` member function, unless
     * preprocessor variable `RANDOMIZE_AS_WORDS` is defined. In that case,
     * A random factor is polled instead.
     *
     * @exception NonRandomizable Thrown if `U` does not support uniform polling
     * and `RANDOMIZE_AS_WORDS` is not defined.
     */
    inline void randomize() {

#ifdef RANDOMIZE_AS_WORDS

        std::vector<FactorTemplate> atoms = atoms();

        *this = atoms[rand() % int(atoms.size())];

#else

        underlying.randomize();

#endif
    }

    /**
     * @brief Generates the atoms.
     *
     * This is a wrapper for the matching `U` member function.
     *
     * @return A vector that contains the atoms.
     */
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

/**
 * @brief Modifies two factors so that they form a left-weighted decomposition.
 *
 * Then returns if something was done.
 *
 * @tparam F A class representing canonical factors.
 * @param u Left factor of the decomposition.
 * @param v Right factor of the decomposition.
 * @return If `u` and `v` were modified.
 */
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

/**
 * @brief Modifies two factors so that they form a right-weighted decomposition.
 *
 * Then returns if something was done.
 *
 * @tparam F A class representing canonical factors.
 * @param u Left factor of the decomposition.
 * @param v Right factor of the decomposition.
 * @return If `u` and `v` were modified.
 */
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

/**
 * @brief Prints factor `f` to output stream `os`.
 *
 * This partially specializes `<<` for `FactorTemplate` classes, as syntactic
 * sugar for `FactorTemplate::print()`.
 *
 * @tparam U A template class for the internal representation of factors.
 * @param os The `IndentedOStream` the factor should be printed in.
 * @param f The factor to be printed.
 * @return A reference to `os`, so that `<<` may be chained.
 */
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
    /**
     * @brief Factor iterator.
     */
    using FactorItr = typename std::list<F>::iterator;

    /**
     * @brief Reverse factor iterator.
     */
    using RevFactorItr = typename std::list<F>::reverse_iterator;

    /**
     * @brief Constant factor iterator.
     */
    using ConstFactorItr = typename std::list<F>::const_iterator;

    /**
     * @brief Constant reverse factor iterator.
     */
    using ConstRevFactorItr = typename std::list<F>::const_reverse_iterator;

    /**
     * @brief Iterator to the first factor.
     *
     * Where the first factor is interpreted as the left-most one.
     *
     * @return An iterator to the first factor.
     */
    inline FactorItr begin() { return factor_list.begin(); }

    /**
     * @brief Reverse iterator to the last factor.
     *
     * Where the last factor is interpreted as the right-most one.
     *
     * @return A reverse iterator to the last factor.
     */
    inline RevFactorItr rbegin() { return factor_list.rbegin(); }

    /**
     * @brief Constant iterator to the first factor.
     *
     * Where the first factor is interpreted as the left-most one.
     *
     * @return A constant iterator to the first factor.
     */
    inline ConstFactorItr cbegin() const { return factor_list.begin(); }

    /**
     * @brief Constant reverse iterator to the last factor.
     *
     * Where the last factor is interpreted as the right-most one.
     *
     * @return A constant reverse iterator to the last factor.
     */
    inline ConstRevFactorItr crbegin() const { return factor_list.rbegin(); }

    /**
     * @brief Iterator to the after-last factor.
     *
     * Where the last factor is interpreted as the right-most one.
     *
     * @return An iterator to the after-last factor.
     */
    inline FactorItr end() { return factor_list.end(); }

    /**
     * @brief Reverse iterator to the before-first factor.
     *
     * Where the first factor is interpreted as the left-most one.
     *
     * @return A reverse iterator to the before-first factor.
     */
    inline RevFactorItr rend() { return factor_list.rend(); }

    /**
     * @brief Constant iterator to the after-last factor.
     *
     * Where the last factor is interpreted as the right-most one.
     *
     * @return A constant iterator to the after-last factor.
     */
    inline ConstFactorItr cend() const { return factor_list.end(); }

    /**
     * @brief Constant reverse iterator to the before-first factor.
     *
     * Where the first factor is interpreted as the left-most one.
     *
     * @return A constant reverse iterator to the before-first factor.
     */
    inline ConstRevFactorItr crend() const { return factor_list.rend(); }

    /**
     * @brief Construct a new BraidTemplate, with a group parameter.
     *
     * @param parameter Group parameter.
     */
    BraidTemplate(Parameter parameter)
        : parameter(parameter), delta(0), factor_list() {}

    /**
     * @brief Construct a new BraidTemplate, from a factor.
     *
     * @param f Factor to be converted to a braid.
     */
    BraidTemplate(const F &f)
        : parameter(f.get_parameter()), delta(0), factor_list() {
        if (f.is_delta()) {
            delta = 1;
        } else if (!f.is_identity()) {
            factor_list.push_back(f);
        }
    }

    /**
     * @brief Converts a string to a `Parameter`.
     *
     * This is a wrapper for the matching `F` static member function.
     *
     * @param str The string to be converted.
     * @return A parameter matching `str`.
     * @exception InvalidStringError Raised if `str` cannot be converted to a
     * `Parameter`.
     */
    inline static Parameter parameter_of_string(const std::string &str) {
        return F::parameter_of_string(str);
    }

    /**
     * @brief Returns the `Parameter` of the braid.
     *
     * @return The `Parameter` of the braid.
     */
    inline Parameter get_parameter() const { return parameter; }

    /**
     * @brief Left-multiplies by the `(delta - inf())`-th power of the Garside
     * element.
     *
     * This is done internally by changing the `delta` member. Therefore, if the
     * braid is in RCF, it corresponds to right multiplication instead.
     *
     * At any rates, if the braid is in a canonical form, `delta` will be its
     * new infimum.
     *
     * @param delta The new infimum.
     */
    inline void set_delta(sint16 delta) { (*this).delta = delta; }

    /**
     * @brief Prints `*this` to `os`.
     *
     * Prints `*this` to `os`, in a format that is compatible with
     * `of_string()` (assuming that the same holds for `F`).
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

    /**
     * @brief Sets the braid to the identity.
     */
    inline void identity() {
        delta = 0;
        factor_list.clear();
    }

    /**
     * @brief Returns the canonical length.
     *
     * _I.e._ the number of factors.
     *
     * @return The canonical length of `*this`.
     */
    inline size_t canonical_length() const { return factor_list.size(); }

    /**
     * @brief Returns the infimum.
     *
     * _I.e._ the exponent of the greatest power of the Garside element that
     * divides the braid.
     *
     * @return The infimum of `*this`.
     */
    inline sint32 inf() const { return delta; }

    /**
     * @brief Returns the supremum.
     *
     * _I.e._ the exponent of the smallest power of the Garside element that
     * the braid divides.
     *
     *  @return The supremum of `*this`.
     */
    inline sint32 sup() const { return inf() + int(canonical_length()); }

    /**
     * @brief Equality check.
     *
     * Assumes that both operands are in the same canonical form.
     *
     * @param v Second operand.
     * @return if `*this` and `v` are equal.
     */
    inline bool compare(const BraidTemplate &v) const {
        return (delta == v.delta && factor_list == v.factor_list);
    }

    /**
     * @brief Equality check.
     *
     * Assumes that both operands are in the same canonical form.
     *
     * Syntactic sugar for `compare()`.
     *
     * @param v Second operand.
     * @return if `*this` and `v` are equal.
     */
    bool operator==(const BraidTemplate &v) const { return compare(v); }

    /**
     * @brief Unequality check.
     *
     * Assumes that both operands are in the same canonical form.
     *
     * @param v Second operand.
     * @return if `*this` and `v` are not equal.
     */
    bool operator!=(const BraidTemplate &v) const { return !compare(v); }

    /**
     * @brief Triviality check.
     *
     * Assumes that the braid is in a canonical form.
     *
     * @return If `*this` is the identity.
     */
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

    /**
     * @brief Computes left meets.
     *
     * @param v Second operand.
     * @return The left meet of `*this` and `v`.
     */
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

    /**
     * @brief Computes left meets, with factors.
     *
     * @param f Second operand (a factor).
     * @return The left meet of `*this` and `f`.
     */
    inline BraidTemplate left_meet(const F &f) {
        return left_meet(BraidTemplate(f));
    }

    /**
     * @brief Computes left meets.
     *
     * Syntactic sugar for `left_meet()`.
     *
     * @param v Second operand.
     * @return The left meet of `*this` and `v`.
     */
    inline BraidTemplate operator^(const BraidTemplate &v) {
        return left_meet(v);
    }

    /**
     * @brief Computes left meets, with factors.
     *
     * Syntactic sugar for `left_meet()`. Notice that the the factor should be the
     * right operand.
     *
     * @param f Second operand (a factor).
     * @return The left meet of `*this` and `f`.
     */
    inline BraidTemplate operator^(const F &f) { return left_meet(f); }

    /**
     * @brief Computes left joins.
     *
     * @param v Second operand.
     * @return The left join of `*this` and `v`.
     */
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

    /**
     * @brief Computes left joins, with factors.
     *
     * @param f Second operand (a factor).
     * @return The left join of `*this` and `f`.
     */
    inline BraidTemplate left_join(const F &f) const {
        return left_join(BraidTemplate(f));
    }

    /**
     * @brief Computes right meets.
     *
     * @param v Second operand.
     * @return The right meet of `*this` and `v`.
     */
    inline BraidTemplate right_meet(const BraidTemplate &v) const {
        return !((!(*this)).left_join(!v));
    }

    /**
     * @brief Computes right meets, with factors.
     *
     * @param f Second operand (a factor).
     * @return The right meet of `*this` and `f`.
     */
    inline BraidTemplate right_meet(const F &f) const {
        return !((!(*this)).left_join(!BraidTemplate(f)));
    }

     /**
     * @brief Computes right joins.
     *
     * @param v Second operand.
     * @return The right join of `*this` and `v`.
     */
    inline BraidTemplate right_join(const BraidTemplate &v) const {
        return !((!(*this)).left_meet(!v));
    }

    /**
     * @brief Computes right joins, with factors.
     *
     * @param f Second operand (a factor).
     * @return The right join of `*this` and `f`.
     */
    inline BraidTemplate right_join(const F &f) const {
        return !((!(*this)).left_meet(!BraidTemplate(f)));
    }

    /**
     * @brief Conjugates by a factor.
     * 
     * @param f The conjugating factor. 
     */
    inline void conjugate(const F &f) {
        left_divide(f);
        right_multiply(f);
    }

    /**
     * @brief Conjugates by a braid.
     * 
     * @param v The conjugating braid. 
     */
    inline void conjugate(const BraidTemplate &v) {
        left_divide(v);
        right_multiply(v);
    }

    /**
     * @brief Conjugates by a factor in RCF.
     * 
     * @param f The conjugating factor. 
     */
    inline void conjugate_rcf(const F &f) {
        left_divide_rcf(f);
        right_multiply_rcf(f);
    }

    /**
     * @brief Conjugates by a braid in RCF.
     * 
     * @param v The conjugating braid. 
     */
    inline void conjugate_rcf(const BraidTemplate &v) {
        left_divide_rcf(v);
        right_multiply_rcf(v);
    }

    /**
     * @brief Returns the first (non-\f$\Delta\f$) factor.
     *
     * Returns the first (non-\f$\Delta\f$) factor. If canonical length is zero,
     * returns the identity factor instead.
     *
     * @return The first (non-\f$\Delta\f$) factor.
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
     * @return The initial factor.
     */
    inline F initial() const { return first().delta_conjugate(-inf()); }

    /**
     * @brief Returns the final factor.
     *
     * Returns the final factor, (i.e. the last one in `factor_list`). If
     * canonical length is zero, returns the identity factor instead.
     *
     * @return The final factor.
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
     * @exception InvalidStringError: Thrown when it isn't possible to
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

/**
 * @brief Prints braid `b` to output stream `os`.
 *
 * This partially specializes `<<` for `BraidTemplate` classes, as syntactic
 * sugar for `BraidTemplate::print()`.
 *
 * @tparam F A template class representing factors.
 * @param os The `IndentedOStream` the braid should be printed in.
 * @param b The braid to be printed.
 * @return A reference to `os`, so that `<<` may be chained.
 */
template <class F>
inline IndentedOStream &operator<<(IndentedOStream &os,
                                   const BraidTemplate<F> &b) {
    b.print(os);
    return os;
}

} // namespace garcide

/**
 * @brief Hash `struct` for `garcide::FactorTemplate`.
 * 
 * Partially specializes `std::hash` for a `garcide::FactorTemplate`.
 * 
 * @tparam U A template class representing the internal data structure of
 * factors.
 */
template <class U> struct std::hash<garcide::FactorTemplate<U>> {
    /**
     * @brief Hash operator.
     * 
     * @param f The factor to be hashed.
     * @return The hash.
     */
    std::size_t operator()(garcide::FactorTemplate<U> const &f) const noexcept {
        return f.hash();
    }
};

/**
 * @brief Hash `struct` for `garcide::BraidTemplate`.
 * 
 * Partially specializes `std::hash` for a `garcide::BraidTemplate`.
 * 
 * @tparam F A template class representing factors.
 */
template <class F> struct std::hash<garcide::BraidTemplate<F>> {
    /**
     * @brief Hash operator.
     * 
     * @param u the braid to be hashed.
     * @return The hash.
     */
    std::size_t operator()(garcide::BraidTemplate<F> const &u) const noexcept {
        return u.hash();
    }
};

#endif
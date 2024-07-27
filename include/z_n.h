#ifndef Z_N
#define Z_N

#include "cgarside.h"

namespace cgarside::z_n {
/**
 * @brief A class for Z ^ n canonical factors.
 *
 * A class for Z ^ n canonical factors (elements of (Z / 2Z) ^ n).
 */
class Underlying {

  public:
    typedef size_t ParameterType;

  private:
    /**
     * @brief dimension.
     *
     * dimension. The n in Z ^ n, represented as a `size_t` (unsigned integer
     * type used for sizes).
     */
    ParameterType dimension;

    /**
     * @brief The factor's coordinates.
     *
     * The factor's coordinates, in Z ^ n canonical basis.
     * As these are all 0 or 1, we use booleans, as (depending on
     * implementation) it should be more space efficient (one bit vers one
     * byte).
     */
    std::vector<bool> coordinates;

  public:
    static ParameterType parameter_of_string(const std::string &str);

    /**
     * @brief Get dimension.
     *
     * Get `dimension` member (which is `private`).
     *
     * @return `dimension`.
     */
    ParameterType GetParameter() const;

    /**
     * @brief Access i-th coordinate.
     *
     * Access i-th coordinate (so that non member functions may read
     * coordinates).
     *
     * @param i The index that is being accessed.
     * @return sint16
     */
    inline bool at(size_t i) const { return coordinates[i]; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `n` as
     * `dimension`.
     *
     * `coordinates` will have length `n`, and will be filled with `false`.
     *
     * @param n The factor's `dimension`.
     */
    Underlying(ParameterType n);

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, tries to extract an
     * atom to set `this` to. If it succeeds, increases `pos` so that it points
     * to just after what was extracted.
     *
     * Letting `Z = '-'? (['1' - '9'] ['0' - '9']* | '0')` be the language of
     * integers, accepted strings are those represented by regular expression
     * `('e' '_'?)? Z`, under the additional hypothesis that the integer they
     * represent is in [`0`, `dimension`[, ignoring whitespaces.
     *
     * "e_i" (with optional "e_") stands for base vector e_i.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception `InvalidStringError`: Thrown when there is no subword starting
     * from `pos` that matches `Z`, or if there is one, if the corresponding
     * integer does not belong to [`1`, `Parameter`[.
     */
    void OfString(const std::string &str, size_t &pos);

    /**
     * @brief Height of the lattice.
     *
     * Height (i.e. `Delta`'s length as a word in the generators, here
     * `dimension`).
     *
     * @return sint16
     */
    sint16 LatticeHeight() const;

    /**
     * @brief Prints internal data in `os`.
     *
     * Prints the factor's `dimension` and `coordinates` to `os`, typically for
     * debugging purposes.
     *
     * @param os The output stream it prints to.
     */
    void Debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor in `os`.
     *
     * Prints the factor in `os` as a product (multiplicative convention) of
     * base vectors.
     *
     * @param os The output stream it prints to.
     */
    void Print(IndentedOStream &os) const;

    /**
     * @brief Sets the factor to the identity.
     *
     * Sets the factor to the identity (i.e. sets all coordinates to `false`).
     */
    void Identity();

    /**
     * @brief Sets the factor to delta.
     *
     * Sets the factor to delta (i.e. sets all coordinates to `true`).
     */
    void Delta();

    /**
     * @brief Computes meets.
     *
     * Computes the (left, although that does not matter here as Z ^ n is
     * abelian) meet of `*this` and `b` (i.e., coordinates-wise `&&`).
     *
     * @param b Second argument.
     * @return The meet of `*this` and `b`.
     */
    Underlying LeftMeet(const Underlying &b) const;

    /**
     * @brief Computes meets.
     *
     * Computes the (right, although that does not matter here as Z ^ n is
     * abelian) meet of `*this` and `b` (i.e., coordinates-wise `&&`).
     *
     * @param b Second argument.
     * @return The meet of `*this` and `b`.
     */
    inline Underlying RightMeet(const Underlying &b) const {
        return LeftMeet(b);
    };

    /**
     * @brief Equality check.
     *
     * Compares `*this` and `b`, returning `true` if they are equal (i.e.
     * coordinate-wise equal).
     *
     * @param b Second argument.
     * @return If `*this` and `b` are equal.
     */
    bool Compare(const Underlying &b) const;

    /**
     * @brief Product computations.
     *
     * Computes the product of `*this` and `b` (i.e. coordinate-wise xor), under
     * the assumption that it lies below delta (this is not actually checked).
     *
     * @param b Second argument.
     * @return The product of `*this` and `b`.
     */
    Underlying Product(const Underlying &b) const;

    /**
     * @brief Complement computations.
     *
     * Computes the complement of `*this` to `b` (i.e. the factor `c` such
     * that `ac=b`), under the assumption that `*this` is smaller than `b`. As Z
     * ^ n is abelian, left and right variants are in fact the same.
     *
     * In this case, this is actually the same (well, as an operation on bit
     * vectors) as product.
     *
     * @param b Second argument.
     * @return The complement of `*this` to `b`.
     */
    inline Underlying LeftComplement(const Underlying &b) const {
        return Product(b);
    };

    /**
     * @brief Complement computations.
     *
     * Computes the complement of `*this` to `b` (i.e. the factor `c` such
     * that `ac=b`), under the assumption that `*this` is smaller than `b`. As Z
     * ^ n is abelian, left and right variants are in fact the same.
     *
     * In this case, this is actually the same (well, as an operation on bit
     * vectors) as product.
     *
     * @param b Second argument.
     * @return The complement of `*this` to `b`.
     */
    inline Underlying RightComplement(const Underlying &b) const {
        return Product(b);
    };

    /**
     * @brief Sets `*this` to a random factor.
     *
     * Sets `*this` to a random factor, following an uniform distribution other
     * factors.
     */
    void Randomize();

    /**
     * @brief List of the atoms.
     *
     * Returns the list of the atoms (i.e. base vectors).
     *
     * @return A vector containing the atoms.
     */
    std::vector<Underlying> Atoms() const;

    /**
     * @brief Conjugates by Delta ^ k.
     *
     * Conjugates `*this` by Delta ^ k (actually, does nothing as Z ^ n is
     * abelian).
     *
     * @param k The exponent.
     */
    inline void DeltaConjugate(sint16 k) {};

    /**
     * @brief Hashes the factor.
     * 
     * Hashes the factor. Done by interpreting it as a polynomial and evaluating it in 2 (yields a bijection between the set of factors, and [0, 2 ^ n[).
     *
     * @return The hash.
     */
    size_t Hash() const;
};

/**
 * @brief Class for canonical factors.
 * 
 * Class for canonical factors, will all corresponding methods.
 */
typedef FactorTemplate<Underlying> Factor;

/**
 * @brief Class for elements of Z ^ n.
 * 
 * Class for elements of Z ^ n, as a Garside group.
 */
typedef BraidTemplate<Factor> Braid;

} // namespace cgarside::z_n

#endif
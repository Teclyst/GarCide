#ifndef ARTIN
#define ARTIN

#include "cgarside.h"
#include "garsideuss.h"

namespace cgarside {

namespace artin {

/// A class for the underlying objects for canonical factors
/// in the Artin presentation braid group case.
/// In this case, permutations.
class Underlying {

  public:
    using Parameter = sint16;

  private:
    Parameter number_of_strands;

    std::vector<sint16> permutation_table;

    /**
     * @brief Maximum braid index.
     *
     * The greatest index that may be used for braids.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const sint16 MAX_NUMBER_OF_STRANDS = 256;

  public:

    static Parameter parameter_of_string(const std::string &str);

    /**
     * @brief Get the `number_of_strands`.
     *
     * Get the `number_of_strands` (which is `private`).
     *
     * @return `number_of_strands`.
     */
    Parameter get_parameter() const;

    sint16 at(size_t i) const { return permutation_table[i]; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `n` as its
     * `number_of_strands`.
     *
     * Its `permutation_table` will have length `n`, and will be filled with
     * zeros (thus this is not a valid factor). Initialize it with `identity`,
     * `delta`, or another similar method.
     *
     * @param n The `number_of_strands` of the factor (also the length of
     * its `permutation_table`).
     */
    Underlying(Parameter n);

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting `Z = -? ([1 - 9] [0 - 9]* | 0)` be the language of integers,
     * accepted strings are those represented by regular expression `Z`, under
     * the additional hypothesis that the integer they represent is in [`1`,
     * `Parameter`[.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception `InvalidStringError`: Thrown when there is no subword starting
     * from `pos` that matches `Z`, or if there is one, if the corresponding
     * integer does not belong to [`1`, `Parameter`[.
     */
    void of_string(const std::string &str, size_t &pos);

    sint16 lattice_height() const;

    /**
     * @brief Prints internal representation to `os`.
     *
     * Prints the factor's `permutation_table` to `os`, typically for debugging
     * purposes.
     *
     * @param os The output stream it prints to.
     */
    void debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints the factor to `os` as a product of atoms.
     *
     * @param os The output stream it prints to.
     */
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    // Compute the meet r of the two factors a and b.  A factor is
    // given as the associated permutation, which is viewed as a
    // bijection on the set {1,...n} and represented as an array whose
    // i-th entry is the image of i under the inverse of the
    // permutation (this convention is different from that in the
    // AsiaCrypt 2001 paper of the author).  The range of indices is
    // [1,n], not [0,n[.  We use a C style array of size (n+1) to
    // represent an n-permutation (the first entry is not used).

    // We define the left meet of two factors a and b to be the
    // longest factor r such that a=ra' and b=rb' for some factors a'
    // and b'.  This coincides with the convention of the paper of
    // Birman, Ko, and Lee, but different from that of the article of
    // Thurston (in Epstein's book).  Indeed, Thurston's is the
    // "right" meet in our sense.
    Underlying left_meet(const Underlying &b) const;

    Underlying right_meet(const Underlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool compare(const Underlying &b) const;

    // product under the hypothesis that it is still simple.
    Underlying product(const Underlying &b) const;

    // Under the assumption a <= b, a.left_complement(b) computes
    // The factor c such that ac = b.
    Underlying left_complement(const Underlying &b) const;

    Underlying right_complement(const Underlying &b) const;

    // Generate a random factor.
    void randomize();

    // List of atoms.
    std::vector<Underlying> atoms() const;

    // Conjugate by delta^k.
    // Used to speed up calculations compared to the default implementation.
    void delta_conjugate_mut(sint16 k);

    size_t hash() const;

    /**
     * @brief Computes the tableau associated with a factor.
     *
     * Computes the tableau associated with `this` and store it in `tab`.
     *
     * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses'
     * code.
     *
     * @param tab A matrix where the tableau is to be stored.
     */
    void tableau(sint16 **&tab) const;

  private:
    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    Underlying inverse() const;

    // Subroutine called by left_meet() and right_meet().
    static void MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s,
                        sint16 t);
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

/**
 * @brief Enum for Thurston types.
 *
 * An enum whose elements represent the three Thurston types.
 */
enum class ThurstonType { Periodic, Reducible, PseudoAsonov };

/**
 * @brief Determines if a braid preserves a family of circles.
 *
 * Determines if braid `b` preserves a family of circles.
 *
 * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses' code.
 *
 * @param b The braid to be tested.
 * @return Whether `b` preserves a family of circles.
 */
bool preserves_circles(const Braid &b);

/**
 * @brief Computes the Thurston type of a braid whose USS was already computed.
 *
 * Computes the Thurston type of braid `b`, using its USS `uss`, in the case
 * where we have access to it.
 *
 * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses' code.
 *
 * @param b The braid whose Thurston type is to be computed.
 * @param uss `b`'s USS.
 * @return `b`'s Thurston type.
 */
ThurstonType thurston_type(const Braid &b,
                           const ultra_summit::UltraSummitSet<Braid> &uss);

/**
 * @brief Computes the Thurston type of a braid.
 *
 * Computes the Thurston type of braid `b`.
 *
 * This was directly copied (mutatis mutandis) from Juan Gonzalez-Meneses' code.
 *
 * @param b The braid whose Thurston type is to be computed.
 * @return `b`'s Thurston type.
 */
ThurstonType thurston_type(const Braid &b);

} // namespace artin

/**
 * @brief Inserts a `ThurstonType` value in the stream.
 *
 * Inserts an `ThurstonType` value in the stream, with obvious conversion
 * conventions.
 *
 * @param type A ThurstonType value.
 * @return A reference to `*this`, so that `<<` may be chained.
 */
template <>
IndentedOStream &IndentedOStream::operator<< <artin::ThurstonType>(
    const artin::ThurstonType &type);

} // namespace cgarside

#endif
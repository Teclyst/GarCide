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

  private:
    sint16 PresentationParameter;

    std::vector<sint16> PermutationTable;

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
    static const sint16 MaxBraidIndex = 256;

  public:
    typedef sint16 ParameterType;

    static ParameterType parameter_of_string(const std::string &str);

    /**
     * @brief Get the `PresentationParameter`.
     *
     * Get the `PresentationParameter` (which is `private`).
     *
     * @return `PresentationParameter`.
     */
    ParameterType GetParameter() const;

    sint16 At(sint16 i) const { return PermutationTable[i]; }
    sint16 &At(sint16 i) { return PermutationTable[i]; }

    /**
     * @brief Construct a new `Underlying`.
     *
     * Construct a new `Underlying`, with `n` as its
     * `PresentationParameter`.
     *
     * Its `PermutationTable` will have length `n`, and will be filled with
     * zeros (thus this is not a valid factor). Initialize it with `Identity`,
     * `Delta`, or another similar method.
     *
     * @param n The `PresentationParameter` of the factor (also the length of
     * its `PermutationTable`).
     */
    Underlying(ParameterType n);

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
    void OfString(const std::string &str, size_t &pos);

    sint16 LatticeHeight() const;

    /**
     * @brief Prints internal representation to `os`.
     *
     * Prints the factor's `PermutationTable` to `os`, typically for debugging
     * purposes.
     *
     * @param os The output stream it prints to.
     */
    void Debug(IndentedOStream &os) const;

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints the factor to `os` as a product of atoms.
     *
     * @param os The output stream it prints to.
     */
    void Print(IndentedOStream &os) const;

    // Set to the Identity element (here the identity).
    void Identity();

    // Set to delta.
    void Delta();

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
    Underlying LeftMeet(const Underlying &b) const;

    Underlying RightMeet(const Underlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool Compare(const Underlying &b) const;

    // Product under the hypothesis that it is still simple.
    Underlying Product(const Underlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    Underlying LeftComplement(const Underlying &b) const;

    Underlying RightComplement(const Underlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<Underlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    void DeltaConjugate(sint16 k);

    size_t Hash() const;

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
    Underlying Inverse() const;

    // Subroutine called by LeftMeet() and RightMeet().
    static void MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s,
                        sint16 t);
};

typedef Factor<Underlying> ArtinFactor;

typedef Braid<ArtinFactor> ArtinBraid;

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
bool preserves_circles(const ArtinBraid &b);

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
ThurstonType thurston_type(const ArtinBraid &b,
                           const ultra_summit::UltraSummitSet<ArtinBraid> &uss);

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
ThurstonType thurston_type(const ArtinBraid &b);

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
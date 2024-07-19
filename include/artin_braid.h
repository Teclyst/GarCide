#include "cgarside.h"

namespace CGarside {

/// A class for the underlying objects for canonical factors
/// in the Artin presentation braid group case.
/// In this case, permutations.
class ArtinBraidUnderlying {

  protected:
    sint16 PresentationParameter;

    std::vector<sint16> PermutationTable;

  public:
    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    sint16 At(sint16 i) const { return PermutationTable[i]; }
    sint16 &At(sint16 i) { return PermutationTable[i]; }

    // Constructor
    ArtinBraidUnderlying(ParameterType n);

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
    ArtinBraidUnderlying LeftMeet(const ArtinBraidUnderlying &b) const;

    ArtinBraidUnderlying RightMeet(const ArtinBraidUnderlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool Compare(const ArtinBraidUnderlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    ArtinBraidUnderlying Inverse() const;

    // Product under the hypothesis that it is still simple.
    ArtinBraidUnderlying Product(const ArtinBraidUnderlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    ArtinBraidUnderlying LeftComplement(const ArtinBraidUnderlying &b) const;

    ArtinBraidUnderlying RightComplement(const ArtinBraidUnderlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<ArtinBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    void DeltaConjugate(sint16 k);

    std::size_t Hash() const {
        std::size_t h = 0;
        for (CGarside::sint16 i = 1; i <= GetParameter(); i++) {
            h = h * 31 + PermutationTable[i];
        }
        return h;
    }

  private:
    // Subroutine called by LeftMeet() and RightMeet().
    static void MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s,
                        sint16 t);
};

typedef Factor<ArtinBraidUnderlying> ArtinBraidFactor;

typedef Braid<ArtinBraidFactor> ArtinBraid;

} // namespace CGarside
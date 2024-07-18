#include "cgarside.h"

namespace CGarside {

class BandBraidUnderlying {

  protected:
    sint16 PresentationParameter;

    sint16 *PermutationTable;

  public:
    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    BandBraidUnderlying(sint16 n);

    BandBraidUnderlying(const BandBraidUnderlying &a);

    ~BandBraidUnderlying() { delete[] PermutationTable; }

    /**
     * @brief Extraction from string.
     *
     * Reads the string `str`, starting at position `pos`, and sets `this` to
     * the corresponding atom.
     *
     * Letting `W = (\s | \t)*` be the language of whitespaces and `Z = -? ([1 -
     * 9] [0 - 9]* | 0)` be the language of integers, accepted strings are those
     * represented by regular expression `\(W Z W,? W Z W\)`, under the
     * additional hypothesis that the two integers lie in [`1`, `Parameter`] and
     * are not equal.
     *
     * @param str The string to extract from.
     * @param pos The position to start from.
     * @exception `InvalidStringError`: Thrown when there is no subword starting
     * from `pos` that matches `\(W Z W,? W Z W\)`, or if there is one, if
     * either integer does not belong to [`1`, `Parameter`], or both are equal.
     */
    void OfString(const std::string &str, size_t &pos);

    /**
     * @brief Prints internal representation to `os`.
     *
     * Prints the factor's `PermutationTable` to `os`, typically for debugging
     * purposes.
     *
     * @param os The output stream it prints to.
     */
    void Debug(IndentedOStream &os) const;

    void AssignDCDT(sint16 *x) const;

    void OfDCDT(const sint16 *x);

    BandBraidUnderlying &Assign(const BandBraidUnderlying &a);

    BandBraidUnderlying &operator=(const BandBraidUnderlying &a);

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints the factor to `os` as a product of atoms.
     *
     * @param os The output stream it prints to.
     */
    void Print(std::ostream &os) const;

    // Set to the Identity element (here the identity).
    void Identity();

    // Set to delta.
    void Delta();

    BandBraidUnderlying LeftMeet(const BandBraidUnderlying &b) const;

    BandBraidUnderlying RightMeet(const BandBraidUnderlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool Compare(const BandBraidUnderlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    BandBraidUnderlying Inverse() const;

    // Product under the hypothesis that it is still simple.
    BandBraidUnderlying Product(const BandBraidUnderlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    BandBraidUnderlying LeftComplement(const BandBraidUnderlying &b) const;

    BandBraidUnderlying RightComplement(const BandBraidUnderlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<BandBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    BandBraidUnderlying DeltaConjugate(sint16 k) const;

    std::size_t Hash() const {
        std::size_t h = 0;
        for (CGarside::sint16 i = 1; i <= GetParameter(); i++) {
            h = h * 31 + PermutationTable[i];
        }
        return h;
    }
};

typedef Factor<BandBraidUnderlying> BandBraidFactor;

typedef Braid<BandBraidFactor> BandBraid;

} // namespace CGarside
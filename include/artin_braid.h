#include "cgarside.h"

namespace CGarside
{

  // A class for the underlying objects for canonical factors
  // in the Artin presentation braid group case.
  // In this case, permutations.
  class ArtinBraidUnderlying
  {

  protected:
    sint16 PresentationParameter;

    sint16 *PermutationTable;

  public:

    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    // Constructor
    ArtinBraidUnderlying(ParameterType n);

    ArtinBraidUnderlying(const ArtinBraidUnderlying &a);

    ~ArtinBraidUnderlying()
    {
      delete[] PermutationTable;
    }

    void OfString(const std::string &str);

    sint16 LatticeHeight() const;

    ArtinBraidUnderlying &Assign(const ArtinBraidUnderlying &a);

    ArtinBraidUnderlying &operator=(const ArtinBraidUnderlying &a);

    void Debug(std::ostream &os) const {
      os << "[";
      for (sint16 i = 1; i <= GetParameter(); i++) {
        os << PermutationTable[i] << " ";
      }
      os << "]";
    }

    // Print to os. Be wary, as it side-effects!
    void Print(std::ostream &os) const;

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
    ArtinBraidUnderlying DeltaConjugate(sint16 k) const;

  private:
    // Subroutine called by LeftMeet() and RightMeet().
    static void MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s, sint16 t);
  };

  typedef Factor<ArtinBraidUnderlying> ArtinBraidFactor;

  typedef Braid<ArtinBraidFactor> ArtinBraid;

}
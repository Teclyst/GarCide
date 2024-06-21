#include "cgarside.h"

namespace CGarside
{

  class BandBraidUnderlying
  {

  protected:
    sint16 PresentationParameter;

    sint16 *PermutationTable;

  public:
    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    // Constructor
    BandBraidUnderlying(sint16 n);

    BandBraidUnderlying(const BandBraidUnderlying &a);

    ~BandBraidUnderlying()
    {
      delete[] PermutationTable;
    }

    void OfString(std::string &str);

    void Debug(std::ostream &os) const
    {
      os << "[";
      for (sint16 i = 1; i <= GetParameter(); i++)
      {
        os << PermutationTable[i] << " ";
      }
      os << "]";
    }

    void AssignDCDT(sint16 *x) const;

    void OfDCDT(const sint16 *x);

    BandBraidUnderlying &Assign(const BandBraidUnderlying &a);

    BandBraidUnderlying &operator=(const BandBraidUnderlying &a);

    // Print to os. Be wary, as it side-effects!
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
    std::list<BandBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    BandBraidUnderlying DeltaConjugate(sint16 k) const;
  };

  typedef Factor<BandBraidUnderlying> BandBraidFactor;

  typedef Braid<BandBraidFactor> BandBraid;

}
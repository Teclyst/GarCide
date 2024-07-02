#include "cgarside.h"

namespace CGarside
{

  class BDualBraidUnderlying
  {

  protected:
    sint16 PresentationParameter;

    sint16 *PermutationTable;

  public:
    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    BDualBraidUnderlying(sint16 n);

    BDualBraidUnderlying(const BDualBraidUnderlying &a);

    ~BDualBraidUnderlying()
    {
      delete[] PermutationTable;
    }

    void OfString(std::string &str);

    void Debug(std::ostream &os) const
    {
      os << "[";
      for (sint16 i = 1; i <= 2 * GetParameter(); i++)
      {
        os << PermutationTable[i] << " ";
      }
      os << "]";
    }

    void AssignDCDT(sint16 *x) const;

    void OfDCDT(const sint16 *x);

    BDualBraidUnderlying &Assign(const BDualBraidUnderlying &a);

    BDualBraidUnderlying &operator=(const BDualBraidUnderlying &a);

    // Print to os. Be wary, as it side-effects!
    void Print(std::ostream &os) const;

    // Set to the Identity element (here the identity).
    void Identity();

    // Set to delta.
    void Delta();

    BDualBraidUnderlying LeftMeet(const BDualBraidUnderlying &b) const;

    BDualBraidUnderlying RightMeet(const BDualBraidUnderlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool Compare(const BDualBraidUnderlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    BDualBraidUnderlying Inverse() const;

    // Product under the hypothesis that it is still simple.
    BDualBraidUnderlying Product(const BDualBraidUnderlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    BDualBraidUnderlying LeftComplement(const BDualBraidUnderlying &b) const;

    BDualBraidUnderlying RightComplement(const BDualBraidUnderlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<BDualBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    BDualBraidUnderlying DeltaConjugate(sint16 k) const;

    std::size_t Hash() const {
      std::size_t h = 0;
      for (CGarside::sint16 i = 1; i <= 2 * GetParameter(); i++)
      {
        h = h * 31 + PermutationTable[i];
      }
      return h;
    }
  };

  typedef Factor<BDualBraidUnderlying> BDualBraidFactor;

  typedef Braid<BDualBraidFactor> BDualBraid;

}
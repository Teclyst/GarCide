#include "cgarside.h"

namespace CGarside
{

  class IDualBraidUnderlying
  {

  protected:
    sint16 PresentationParameter;

    // 0 for Identity, 1 for Delta, 2 for a reflection.
    sint16 Type;

    // The reflection Axis.
    sint16 Axis;

  public:
    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    IDualBraidUnderlying(sint16 n);

    IDualBraidUnderlying(const IDualBraidUnderlying &a);

    void OfString(std::string &str);

    void Debug(std::ostream &os) const
    {
      os << "(" << Type << ", " << Axis << ") ";
    }

    IDualBraidUnderlying &Assign(const IDualBraidUnderlying &a);

    IDualBraidUnderlying &operator=(const IDualBraidUnderlying &a);

    // Print to os. Be wary, as it side-effects!
    void Print(std::ostream &os) const;

    // Set to the Identity element (here the identity).
    void Identity();

    // Set to delta.
    void Delta();

    IDualBraidUnderlying LeftMeet(const IDualBraidUnderlying &b) const;

    IDualBraidUnderlying RightMeet(const IDualBraidUnderlying &b) const;

    // Equality check.
    // We check wether the underlying permutation table are (pointwise) equal.
    bool Compare(const IDualBraidUnderlying &b) const;

    // Product under the hypothesis that it is still simple.
    IDualBraidUnderlying Product(const IDualBraidUnderlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    IDualBraidUnderlying LeftComplement(const IDualBraidUnderlying &b) const;

    IDualBraidUnderlying RightComplement(const IDualBraidUnderlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<IDualBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    IDualBraidUnderlying DeltaConjugate(sint16 k) const;

    std::size_t Hash() const {
      std::size_t h = Axis;
      return h;
    }
  };

  typedef Factor<IDualBraidUnderlying> IDualBraidFactor;

  typedef Braid<IDualBraidFactor> IDualBraid;

}
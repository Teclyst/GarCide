#include "cgarside.h"

namespace CGarside {

class IDualBraidUnderlying {

  protected:
    sint16 PresentationParameter;

    // 0 for Identity, 1 for Delta, 2 for a reflection.
    sint16 Type;

    // The point the reflection sends 0 on.
    sint16 Point;

  public:
    typedef sint16 ParameterType;

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    IDualBraidUnderlying(sint16 n);

    void OfString(const std::string &str, size_t &pos);

    void Debug(IndentedOStream &os) const {
        os << "(" << Type << ", " << Point << ") ";
    }

    // Print to os.
    void Print(IndentedOStream &os) const;

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
    // Used to speed up calculations compared to default implementation.
    IDualBraidUnderlying DeltaConjugate(sint16 k) const;

    inline std::size_t Hash() const {
        std::size_t h = Point;
        return h;
    }
};

typedef Factor<IDualBraidUnderlying> IDualBraidFactor;

typedef Braid<IDualBraidFactor> IDualBraid;

} // namespace CGarside
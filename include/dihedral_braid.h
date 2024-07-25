#include "cgarside.h"

namespace cgarside::dihedral {

class Underlying {

  protected:
    sint16 PresentationParameter;

    // 0 for Identity, 1 for Delta, 2 for a reflection.
    sint16 Type;

    // The point the reflection sends 0 on.
    sint16 Point;

  public:
    typedef sint16 ParameterType;

    static ParameterType parameter_of_string(const std::string &str);

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    Underlying(sint16 n);

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
    // Used to speed up calculations compared to default implementation.
    Underlying DeltaConjugate(sint16 k) const;

    inline std::size_t Hash() const {
        std::size_t h = Point;
        return h;
    }
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

} // namespace cgarside::dihedral
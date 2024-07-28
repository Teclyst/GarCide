#include "cgarside.h"

namespace cgarside::dihedral {

class Underlying {

  protected:
    sint16 PresentationParameter;

    // 0 for identity, 1 for delta, 2 for a reflection.
    sint16 Type;

    // The point the reflection sends 0 on.
    sint16 Point;

  public:
    typedef sint16 Parameter;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    sint16 lattice_height() const;

    // Constructor
    Underlying(sint16 n);

    void of_string(const std::string &str, size_t &pos);

    void debug(IndentedOStream &os) const {
        os << "(" << Type << ", " << Point << ") ";
    }

    // print to os.
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

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
    // Used to speed up calculations compared to default implementation.
    Underlying delta_conjugate_mut(sint16 k) const;

    inline std::size_t hash() const {
        std::size_t h = Point;
        return h;
    }
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

} // namespace cgarside::dihedral
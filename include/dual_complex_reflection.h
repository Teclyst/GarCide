#include "cgarside.h"

namespace cgarside {

namespace dual_complex {

// We represent B(e, e, n + 1).
struct Parameter {
    sint16 e;
    sint16 n;

    Parameter(sint16 e, sint16 n) : e(e), n(n) {}

    inline bool compare(const Parameter &p) const {
        return ((e == p.e) && (n == p.n));
    }

    inline bool operator==(const Parameter &p) const {
        return compare(p);
    }

    inline bool operator!=(const Parameter &p) const {
        return !compare(p);
    }

    void print(IndentedOStream &os) const;
};

// You may use Braids with parameter e, n such that e * n <= `MaxE` *
// `MaxN`. Note that `MaxE` IS NOT a strict bound; rather, it is the
// greatest possible e such that it is possible to have braids with parameter n
// as great as `MaxN`. To understand what we are doing, it is advised
// to look at `arXiv:math/0403400v2` (Bessis, Corran, 2004).
class Underlying {

  protected:
    Parameter PresentationParameter;

    // The induced permutation, where 0 is the 0-th coordinates, and then the
    // next n coordinates represent the powers w ^ -i (with i ranging between 0
    // and n - 1) of an ne-th root of unity w. We use the same conventions as
    // before: letting sigma be the induced permutation, then
    // permutation_table[i] is sigma^(-1)(i).
    std::vector<sint16> permutation_table;

    // The multiplicating coefficients. These are e-th roots of unity, with the
    // added condition of their product's being one. They are represented by
    // integer ranging between 0 and e - 1, with i standing for w^(ei) (with w
    // the same root as for the permutation_table).
    std::vector<sint16> CoefficientTable;

  public:
    /**
     * @brief Maximum value for `PresentationParameter.n`.
     *
     * The greatest `PresentationParameter.n` value that may be used for braids.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const sint16 MaxN = 256;

    /**
     * @brief Maximum value for `PresentationParameter.e`.
     *
     * Not exactly. This is not a strict bound on `PresentationParameter.e`. The
     * real condition is `PresentationParameter.e * PresentationParameter.n <=
     * MaxE * MaxN`.
     *
     * It is used because we use `thread_local` objects to avoid some
     * allocations, and their size must be known at compile time.
     *
     * Having too big `thread_local` objects might cause some issue with thread
     * spawning.
     */
    static const sint16 MaxE = 1;

    typedef Parameter Parameter;

    static Parameter parameter_of_string(const std::string &str);

    Parameter get_parameter() const;

    sint16 lattice_height() const;

    // Constructor
    Underlying(Parameter p);

    void of_string(const std::string &str, size_t &pos);

    void debug(IndentedOStream &os) const;

    void AssignPartition(sint16 *x) const;

    void OfPartition(const sint16 *x);

    // print to os. Be wary, as it side-effects!
    void print(IndentedOStream &os) const;

    // Set to the identity element (here the identity).
    void identity();

    // Set to delta.
    void delta();

    Underlying
    left_meet(const Underlying &b) const;

    inline Underlying
    right_meet(const Underlying &b) const {
        return left_meet(b);
    };

    // Equality check.
    // We check whether the underlying permutation table are (pointwise) equal.
    bool compare(const Underlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    Underlying Inverse() const;

    // product under the hypothesis that it is still simple.
    Underlying
    product(const Underlying &b) const;

    // Under the assumption a <= b, a.left_complement(b) computes
    // The factor c such that ac = b.
    Underlying
    left_complement(const Underlying &b) const;

    Underlying
    right_complement(const Underlying &b) const;

    // Generate a random factor.
    void randomize();

    // List of atoms.
    std::vector<Underlying> atoms() const;

    // Conjugate by delta^k.
    // Used to speed up calculations compared to the default implementation.
    void delta_conjugate_mut(sint16 k);

    std::size_t hash() const;
};

typedef FactorTemplate<Underlying> Factor;

typedef BraidTemplate<Factor> Braid;

} // namespace CGarside

template <>
IndentedOStream &IndentedOStream::operator<< <dual_complex::Parameter>(
    const dual_complex::Parameter &p);

}
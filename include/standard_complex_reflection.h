#include "cgarside.h"

namespace CGarside {

// We represent B(e, e, n).
struct ComplexStandardBraidParameter {
    sint16 e;
    sint16 n;

    ComplexStandardBraidParameter(sint16 e, sint16 n) : e(e), n(n) {}

    inline bool Compare(const ComplexStandardBraidParameter &p) const {
        return ((e == p.e) && (n == p.n));
    }

    inline bool operator==(const ComplexStandardBraidParameter &p) const {
        return Compare(p);
    }

    inline bool operator!=(const ComplexStandardBraidParameter &p) const {
        return !Compare(p);
    }

    void Debug(std::ostream &os) const {
        os << "{ e: " << e << ", n: " << n << " }";
    }
};

std::ostream &operator<<(std::ostream &os,
                         const ComplexStandardBraidParameter &p);

// You may use Braids with parameter e, n such that e * n <= `MaxE` *
// `MaxBraidIndex`. Note that `MaxE` IS NOT a strict bound; rather, it is the
// greatest possible e such that it is possible to have braids with parameter n
// as great as `MaxBraidIndex`. To understand what we are doing, it is advised
// to look at `arXiv:math/0403400v2` (Bessis, Corran, 2004).
class ComplexStandardBraidUnderlying {

  protected:
    ComplexStandardBraidParameter PresentationParameter;

    // The induced permutation, where 0 is the 0-th coordinates, and then the
    // next n coordinates represent the powers w ^ -i (with i ranging between 0
    // and n - 1) of an ne-th root of unity w. We use the same conventions as
    // before: letting sigma be the induced permutation, then
    // PermutationTable[i] is sigma^(-1)(i).
    std::vector<sint16> PermutationTable;

    // The multiplicating coefficients. These are e-th roots of unity, with the
    // added condition of their product's being one. They are represented by
    // integer ranging between 0 and e - 1, with i standing for w^(ei) (with w
    // the same root as for the PermutationTable).
    std::vector<sint16> CoefficientTable;

  public:
    typedef ComplexStandardBraidParameter ParameterType;

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    ComplexStandardBraidUnderlying(ComplexStandardBraidParameter p);

    void OfString(std::string str);

    // Prints to `os` the internal representation of the factor.
    // @param os The `std::ostream` we are printing to.
    void Debug(std::ostream &os) const {
        os << "[";
        for (sint16 i = 0; i < GetParameter().n; i++) {
            os << "(" << PermutationTable[i] << ", " << CoefficientTable[i]
               << "), ";
        }
        os << "]";
    }

    void AssignPartition(sint16 *x) const;

    void OfPartition(const sint16 *x);

    /**
     * @brief Prints the factor to `os`.
     *
     * Prints to `os` Neaime's shortest word representative in the
     * Corran-Picantin generators (See George Neaime, Interval Garside
     * Structures for the Complex Braid Groups,
     * [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)))
     *
     * @param os The output stream it prints to.
     */
    void Print(std::ostream &os) const;

    // Set to the Identity element (here the identity).
    void Identity();

    // Set to delta.
    void Delta();

    // Gets the direct permutation table associated with the factor.
    // (i.e., perm[i] is the image of i by the permutation).
    // Used as dividing by atoms requires columns to be easily findable.
    // @param dir_perm The table that is to be filled.
    void Direct(sint16 *dir_perm) const;

    // Checks if s_i left divides `*this`. (See George Neaime, Interval Garside
    // Structures for the Complex Braid Groups, Proposition 3.13,
    // [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor.
    // @param i The integer i for which we check if s_i left divides the factor.
    bool IsSLeftDivisor(const sint16 *dir_perm, sint16 i) const {
        return (PermutationTable[i - 1] > PermutationTable[i - 2])
                   ? (CoefficientTable[i - 1] != 0)
                   : (CoefficientTable[i - 2] == 0);
    };

    // Checks if t_i left divides `*this`. (See George Neaime, Interval Garside
    // Structures for the Complex Braid Groups, Proposition 3.13,
    // [arXiv:1707.06864](https://arxiv.org/abs/1707.06864)).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor.
    // @param i The integer i for which we check if t_i left divides the factor.
    bool IsTLeftDivisor(const sint16 *dir_perm, sint16 i) const {
        return (PermutationTable[1] > PermutationTable[0])
                   ? (CoefficientTable[1] != 0)
                   : (CoefficientTable[0] ==
                      ((i == 0) ? 0 : GetParameter().e - i));
    };

    // Left multiplies by s_i (or divides, which is the same as it has order 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by s_i the factor.
    void SLeftMultiply(sint16 *dir_perm, sint16 i) {
        std::swap(CoefficientTable[i - 1], CoefficientTable[i - 2]);
        std::swap(PermutationTable[i - 1], PermutationTable[i - 2]);
        std::swap(dir_perm[PermutationTable[i - 1]],
                  dir_perm[PermutationTable[i - 2]]);
    };

    // Left multiplies by t_i (or divides, which is the same as it has order 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by t_i the factor.
    void TLeftMultiply(sint16 *dir_perm, sint16 i) {
        std::swap(CoefficientTable[0], CoefficientTable[1]);
        std::swap(PermutationTable[0], PermutationTable[1]);
        CoefficientTable[0] = Rem(CoefficientTable[0] - i, GetParameter().e);
        CoefficientTable[1] = Rem(CoefficientTable[1] + i, GetParameter().e);
        std::swap(dir_perm[PermutationTable[0]], dir_perm[PermutationTable[1]]);
    };

    // Right multiplies by s_i (or divides, which is the same as it has order
    // 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by s_i the factor.
    void SRightMultiply(sint16 *dir_perm, sint16 i) {
        std::swap(PermutationTable[dir_perm[i - 1]],
                  PermutationTable[dir_perm[i - 2]]);
        std::swap(dir_perm[i - 1], dir_perm[i - 2]);
    };

    // Right multiplies by t_i (or divides, which is the same as it has order
    // 2).
    // @param dir_perm A table that holds the direct table of the permutation
    // induced by the factor. Modified by this function so that it remains up to
    // date.
    // @param i The integer i for which we multiply by t_i the factor.
    void TRightMultiply(sint16 *dir_perm, sint16 i) {
        std::swap(PermutationTable[dir_perm[0]], PermutationTable[dir_perm[1]]);
        CoefficientTable[dir_perm[0]] =
            Rem(CoefficientTable[dir_perm[0]] - i, GetParameter().e);
        CoefficientTable[dir_perm[1]] =
            Rem(CoefficientTable[dir_perm[1]] + i, GetParameter().e);
        std::swap(dir_perm[0], dir_perm[1]);
    };

    ComplexStandardBraidUnderlying
    LeftMeet(const ComplexStandardBraidUnderlying &b) const;

    inline ComplexStandardBraidUnderlying
    RightMeet(const ComplexStandardBraidUnderlying &b) const {
        return Inverse().LeftMeet(b.Inverse()).Inverse();
    };

    // Equality check.
    // We check whether the underlying permutation table are (pointwise) equal.
    bool Compare(const ComplexStandardBraidUnderlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    ComplexStandardBraidUnderlying Inverse() const;

    // Product under the hypothesis that it is still simple.
    ComplexStandardBraidUnderlying
    Product(const ComplexStandardBraidUnderlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    ComplexStandardBraidUnderlying
    LeftComplement(const ComplexStandardBraidUnderlying &b) const;

    ComplexStandardBraidUnderlying
    RightComplement(const ComplexStandardBraidUnderlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<ComplexStandardBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    void DeltaConjugate(sint16 k);

    std::size_t Hash() const {
        std::size_t h = 0;
        for (CGarside::sint16 i = 1; i < GetParameter().n; i++) {
            h = h * 31 + PermutationTable[i];
        }
        for (CGarside::sint16 i = 1; i < GetParameter().n; i++) {
            h = h * 31 + CoefficientTable[i];
        }
        return h;
    }
};

std::ostream &operator<<(std::ostream &os,
                         const ComplexStandardBraidUnderlying &u);

typedef Factor<ComplexStandardBraidUnderlying> ComplexStandardBraidFactor;

typedef Braid<ComplexStandardBraidFactor> ComplexStandardBraid;

} // namespace CGarside
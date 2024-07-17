#include "cgarside.h"

namespace CGarside
{
  // The greates possible value for e such that it is possible to have braids with parameter n as great as `MaxBraidIndex`
  // It is possible to use objects with e > MaxE, but then the other parameter will have to be smaller than `MaxBraidIndex` / e.
  const sint16 MaxE = 1;

  // We represent B(e, e, n + 1).
  struct ComplexDualBraidParameter
  {
    sint16 e;
    sint16 n;

    ComplexDualBraidParameter(sint16 e, sint16 n) : e(e), n(n) {}

    inline bool Compare(const ComplexDualBraidParameter &p) const
    {
      return ((e == p.e) && (n == p.n));
    }

    inline bool operator==(const ComplexDualBraidParameter &p) const
    {
      return Compare(p);
    }

    inline bool operator!=(const ComplexDualBraidParameter &p) const
    {
      return !Compare(p);
    }

    void Debug(std::ostream &os) const {
      os << "{ e: " << e << ", n: " << n << " }"; 
    }

    void check_non_crossing(sint16 *x) {
      for (sint16 i = 1; i < e * n; i++) {
        for (sint16 j = i + 1; j <= e * n; j++) {
          for (sint16 k = i + 1; k <= j; k++) {
            for (sint16 l = j + 1; l <= e * n; l++) {
              if ((x[i] == x[j]) && (x[k] == x[l]) && (x[i] != x[k])) {
                exit(1);
              }
            }
          }
        }
      }
    }

  };

  std::ostream& operator<<(std::ostream &os, const ComplexDualBraidParameter &p);

  // You may use Braids with parameter e, n such that e * n <= `MaxE` * `MaxBraidIndex`.
  // Note that `MaxE` IS NOT a strict bound; rather, it is the greatest possible e such that it is possible to have braids with parameter n as great as `MaxBraidIndex`.
  // To understand what we are doing, it is advised to look at `arXiv:math/0403400v2` (Bessis, Corran, 2004).
  class ComplexDualBraidUnderlying
  {

  protected:
    ComplexDualBraidParameter PresentationParameter;

    // The induced permutation, where 0 is the 0-th coordinates, and then the next n coordinates represent the powers w ^ -i (with i ranging between 0 and n - 1) of an ne-th root of unity w.
    // We use the same conventions as before: letting sigma be the induced permutation, then PermutationTable[i] is sigma^(-1)(i).
    std::vector<sint16> PermutationTable;

    // The multiplicating coefficients. These are e-th roots of unity, with the added condition of their product's being one.
    // They are represented by integer ranging between 0 and e - 1, with i standing for w^(ei) (with w the same root as for the PermutationTable).
    std::vector<sint16> CoefficientTable;

  public:
    typedef ComplexDualBraidParameter ParameterType;

    ParameterType GetParameter() const;

    sint16 LatticeHeight() const;

    // Constructor
    ComplexDualBraidUnderlying(ComplexDualBraidParameter p);

    void OfString(std::string str);

    void Debug(std::ostream &os) const
    {
      os << "[";
      for (sint16 i = 0; i <= GetParameter().n; i++)
      {
        os << "(" << PermutationTable[i] << ", " << CoefficientTable[i] << "), ";
      }
      os << "]";
    }

    void AssignPartition(sint16 *x) const;

    void OfPartition(const sint16 *x);

    // Print to os. Be wary, as it side-effects!
    void Print(std::ostream &os) const;

    // Set to the Identity element (here the identity).
    void Identity();

    // Set to delta.
    void Delta();

    ComplexDualBraidUnderlying LeftMeet(const ComplexDualBraidUnderlying &b) const;

    inline ComplexDualBraidUnderlying RightMeet(const ComplexDualBraidUnderlying &b) const
    {
      return LeftMeet(b);
    };

    // Equality check.
    // We check whether the underlying permutation table are (pointwise) equal.
    bool Compare(const ComplexDualBraidUnderlying &b) const;

    // Computes the factor corresponding to the inverse permutation.
    // Used to simplify complement operation.
    ComplexDualBraidUnderlying Inverse() const;

    // Product under the hypothesis that it is still simple.
    ComplexDualBraidUnderlying Product(const ComplexDualBraidUnderlying &b) const;

    // Under the assumption a <= b, a.LeftComplement(b) computes
    // The factor c such that ac = b.
    ComplexDualBraidUnderlying LeftComplement(const ComplexDualBraidUnderlying &b) const;

    ComplexDualBraidUnderlying RightComplement(const ComplexDualBraidUnderlying &b) const;

    // Generate a random factor.
    void Randomize();

    // List of atoms.
    std::vector<ComplexDualBraidUnderlying> Atoms() const;

    // Conjugate by Delta^k.
    // Used to speed up calculations compared to the default implementation.
    void DeltaConjugate(sint16 k);

    std::size_t Hash() const
    {
      std::size_t h = 0;
      for (CGarside::sint16 i = 1; i <= GetParameter().n; i++)
      {
        h = h * 31 + PermutationTable[i];
      }
      for (CGarside::sint16 i = 1; i <= GetParameter().n; i++)
      {
        h = h * 31 + CoefficientTable[i];
      }
      return h;
    }
  };

  std::ostream &operator<<(std::ostream &os, const ComplexDualBraidUnderlying &u);

  typedef Factor<ComplexDualBraidUnderlying> ComplexDualBraidFactor;

  typedef Braid<ComplexDualBraidFactor> ComplexDualBraid;

}
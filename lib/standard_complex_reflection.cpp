#include "standard_complex_reflection.h"
#include <string>
#include <iostream>

namespace CGarside
{
  std::ostream &operator<<(std::ostream &os, const ComplexStandardBraidParameter &p)
  {
    p.Debug(os);
    return os;
  };

  ComplexStandardBraidParameter ComplexStandardBraidUnderlying::GetParameter() const
  {
    return PresentationParameter;
  }

  sint16 ComplexStandardBraidUnderlying::LatticeHeight() const
  {
    sint16 n = GetParameter().n;
    return n * (n + 1);
  }

  ComplexStandardBraidUnderlying::ComplexStandardBraidUnderlying(ComplexStandardBraidParameter p)
      : PresentationParameter(p), PermutationTable(p.n), CoefficientTable(p.n) {}

  void ComplexStandardBraidUnderlying::Print(std::ostream &os) const
  {
    thread_local sint16 dir_perm[MaxBraidIndex];
    ComplexStandardBraidUnderlying copy = *this;
    copy.Direct(dir_perm);
    sint16 n = GetParameter().n, e = GetParameter().e;
    for (sint16 i = 2; i <= n; i++)
    {
      for (sint16 j = i; j > 2; j--)
      {
        if (copy.IsSLeftDivisor(dir_perm, j))
        {
          copy.SLeftMultiply(dir_perm, j);
          os << "s" << j << " ";
        }
        else
        {
          goto end_subword;
        }
      }

      if (copy.IsTLeftDivisor(dir_perm, 1))
      {
        copy.TLeftMultiply(dir_perm, 1);
        os << "t" << 1 << " ";

        if (copy.IsTLeftDivisor(dir_perm, 0))
        {
          copy.TLeftMultiply(dir_perm, 0);
          os << "t" << 0 << " ";

          for (sint16 j = 3; j < i + 1; j++)
          {
            if (copy.IsSLeftDivisor(dir_perm, j))
            {
              copy.SLeftMultiply(dir_perm, j);
              os << "s" << j << " ";
            }
            else
            {
              goto end_subword;
            }
          }
        }
        goto end_subword;
      }

      for (sint16 k = 0; k < e; k++)
      {
        if (copy.IsTLeftDivisor(dir_perm, k))
        {
          copy.TLeftMultiply(dir_perm, k);
          os << "t" << k << " ";
          goto end_subword;
        }
      }
    end_subword:;
    }
  }

  void ComplexStandardBraidUnderlying::OfString(std::string str)
  {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::regex factor{"(s|t)(" + number_regex + ")"};
    std::smatch match;
    if (std::regex_search(str, match, factor, std::regex_constants::match_continuous))
    {
      sint16 i = std::stoi(match[2]);
      if (match[1] == "s")
      {
        if ((3 <= i) && (i <= n))
        {
          // s_i.
          Identity();
          PermutationTable[i - 2] = i - 1;
          PermutationTable[i - 1] = i - 2;
        }
        else
        {
          throw InvalidStringError();
        }
      }
      else
      {
        if ((3 <= i) && (i <= n))
        {
          // t_i.
          i = Rem(i, e);
          Identity();
          PermutationTable[0] = 1;
          PermutationTable[1] = 0;
          CoefficientTable[1] = i;
          CoefficientTable[0] = e - ((i == 0) ? e : i);
        }
      }
    }
    else
    {
      throw InvalidStringError();
    }
  }

  void ComplexStandardBraidUnderlying::Direct(sint16 *dir_perm) const
  {
    for (sint16 i = 0; i < GetParameter().n; i++)
    {
      dir_perm[PermutationTable[i]] = i;
    }
  };

  ComplexStandardBraidUnderlying ComplexStandardBraidUnderlying::LeftMeet(const ComplexStandardBraidUnderlying &b) const
  {
    thread_local sint16 dir_perm_a[MaxBraidIndex], dir_perm_b[MaxBraidIndex], dir_perm_meet[MaxBraidIndex];
    ComplexStandardBraidUnderlying a_copy = *this;
    ComplexStandardBraidUnderlying b_copy = b;
    ComplexStandardBraidUnderlying meet(GetParameter());
    meet.Identity();
    a_copy.Direct(dir_perm_a);
    b_copy.Direct(dir_perm_b);
    meet.Direct(dir_perm_meet);
    sint16 n = GetParameter().n, e = GetParameter().e;
    for (sint16 i = 2; i <= n; i++)
    {
      for (sint16 j = i; j > 2; j--)
      {
        if (a_copy.IsSLeftDivisor(dir_perm_a, j) && b_copy.IsSLeftDivisor(dir_perm_b, j))
        {
          meet.SRightMultiply(dir_perm_meet, j);
          a_copy.SLeftMultiply(dir_perm_a, j);
          b_copy.SLeftMultiply(dir_perm_b, j);
        }
        else
        {
          goto end_subword;
        }
      }

      if (a_copy.IsTLeftDivisor(dir_perm_a, 0) && b_copy.IsTLeftDivisor(dir_perm_b, 0))
      {

        meet.TRightMultiply(dir_perm_meet, 0);
        a_copy.TLeftMultiply(dir_perm_a, 0);
        b_copy.TLeftMultiply(dir_perm_b, 0);

        if (a_copy.IsTLeftDivisor(dir_perm_a, e - 1) && b_copy.IsTLeftDivisor(dir_perm_b, e - 1))
        {
          meet.TRightMultiply(dir_perm_meet, e - 1);
          a_copy.TLeftMultiply(dir_perm_a, e - 1);
          b_copy.TLeftMultiply(dir_perm_b, e - 1);

          for (sint16 j = 3; j < i + 1; j++)
          {
            if (a_copy.IsSLeftDivisor(dir_perm_a, j) && b_copy.IsSLeftDivisor(dir_perm_b, j))
            {
              meet.SRightMultiply(dir_perm_meet, j);
              a_copy.SLeftMultiply(dir_perm_a, j);
              b_copy.SLeftMultiply(dir_perm_b, j);
            }
            else
            {
              goto end_subword;
            }
          }
        }
        goto end_subword;
      }

      for (sint16 k = 1; k < e; k++)
      {
        if (a_copy.IsTLeftDivisor(dir_perm_a, k) && b_copy.IsTLeftDivisor(dir_perm_b, k))
        {
          meet.TRightMultiply(dir_perm_meet, k);
          a_copy.TLeftMultiply(dir_perm_a, k);
          b_copy.TLeftMultiply(dir_perm_b, k);
          goto end_subword;
        }
      }
    end_subword:;
    }
    return meet;
  }

  void ComplexStandardBraidUnderlying::Identity()
  {
    for (sint16 i = 0; i < GetParameter().n; i++)
    {
      PermutationTable[i] = i;
      CoefficientTable[i] = 0;
    }
  }

  void ComplexStandardBraidUnderlying::Delta()
  {
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 1; i < n; i++)
    {
      PermutationTable[i] = i;
      CoefficientTable[i] = 1;
    }
    PermutationTable[0] = 0;
    CoefficientTable[0] = Rem(-n + 1, e);
  }

  bool ComplexStandardBraidUnderlying::Compare(const ComplexStandardBraidUnderlying &b) const
  {
    sint16 i;
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter().n, b.GetParameter().n);
    }
    for (i = 0; i < GetParameter().n; i++)
    {
      if ((PermutationTable[i] != b.PermutationTable[i]) || (CoefficientTable[i] != b.CoefficientTable[i]))
      {
        return false;
      }
    }
    return true;
  };

  ComplexStandardBraidUnderlying ComplexStandardBraidUnderlying::Inverse() const
  {
    ComplexStandardBraidUnderlying f = ComplexStandardBraidUnderlying(GetParameter());
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 0; i < n; i++)
    {
      f.PermutationTable[PermutationTable[i]] = i;
      f.CoefficientTable[PermutationTable[i]] = CoefficientTable[i] == 0 ? 0 : e - CoefficientTable[i];
    }
    return f;
  };

  ComplexStandardBraidUnderlying ComplexStandardBraidUnderlying::Product(const ComplexStandardBraidUnderlying &b) const
  {
    ComplexStandardBraidUnderlying f = ComplexStandardBraidUnderlying(GetParameter());
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter().n, b.GetParameter().n);
    }
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 0; i < n; i++)
    {
      f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
      f.CoefficientTable[i] = Rem(b.CoefficientTable[PermutationTable[i]] + CoefficientTable[i], e);
    }
    return f;
  };

  ComplexStandardBraidUnderlying ComplexStandardBraidUnderlying::LeftComplement(const ComplexStandardBraidUnderlying &b) const
  {
    return b.Product(Inverse());
  };

  ComplexStandardBraidUnderlying ComplexStandardBraidUnderlying::RightComplement(const ComplexStandardBraidUnderlying &b) const
  {
    return Inverse().Product(b);
  };

  void ComplexStandardBraidUnderlying::DeltaConjugate(sint16 k)
  {
    // Delta is diagonal, and acts almost homothetically.
    // Therefore conjugating by some power of it does nothing on most coefficients.
    sint16 n = GetParameter().n, e = GetParameter().e;

    // In this case `*this` commutes with Delta.
    if (PermutationTable[0] == 0)
    {
      return;
    }

    // Otherwise the two non trivial coefficient are 0 and the i such that `PermutationTable[i] == 0`.
    CoefficientTable[0] = Rem(CoefficientTable[0] + k * n, e);
    for (sint16 i = 1; i < n; i++)
    {
      if (PermutationTable[i] == 0)
      {
        CoefficientTable[i] = Rem(CoefficientTable[i] - k * n, e);
        return;
      }
    }
  }

  void ComplexStandardBraidUnderlying::Randomize()
  {
    throw NonRandomizable();
  }

  std::vector<ComplexStandardBraidUnderlying> ComplexStandardBraidUnderlying::Atoms() const
  {
    ComplexStandardBraidParameter p = GetParameter();
    std::vector<ComplexStandardBraidUnderlying> atoms;
    ComplexStandardBraidUnderlying atom = ComplexStandardBraidUnderlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 2; i <= n - 1; i++)
    {
      // s_(i+1).
      atom.Identity();
      atom.PermutationTable[i - 1] = i;
      atom.PermutationTable[i] = i - 1;
      atoms.push_back(atom);
    }
    for (sint16 k = 0; k < e; k++)
    {
      // t_k.
      atom.Identity();
      atom.PermutationTable[0] = 1;
      atom.PermutationTable[1] = 0;
      atom.CoefficientTable[1] = k;
      atom.CoefficientTable[0] = e - ((k == 0) ? e : k);
      atoms.push_back(atom);
    }
    return atoms;
  }

  std::ostream &operator<<(std::ostream &os, const ComplexStandardBraidUnderlying &u)
  {
    u.Print(os);
    return os;
  }
}
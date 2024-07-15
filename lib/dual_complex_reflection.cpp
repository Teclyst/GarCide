#include "dual_complex_reflection.h"
#include <string>
#include <iostream>

namespace CGarside
{
  std::ostream &operator<<(std::ostream &os, const ComplexDualBraidParameter &p)
  {
    p.Debug(os);
    return os;
  };

  ComplexDualBraidParameter ComplexDualBraidUnderlying::GetParameter() const
  {
    return PresentationParameter;
  }

  sint16 ComplexDualBraidUnderlying::LatticeHeight() const
  {
    return GetParameter().n + 1;
  }

  ComplexDualBraidUnderlying::ComplexDualBraidUnderlying(ComplexDualBraidParameter p)
      : PresentationParameter(p), PermutationTable(p.n + 1), CoefficientTable(p.n + 1) {}

  void ComplexDualBraidUnderlying::Print(std::ostream &os) const
  {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::vector<bool> seen(n + 1, false);
    std::vector<sint16> curr_cycle;
    sint16 other_smallest = 0, cycle_type, curr, c = 0;
    seen[0] = true;

    curr = PermutationTable[0];
    // Short assymetric case.
    if (curr != 0)
    {
      while (curr != 0)
      {
        curr_cycle.push_back(curr);
        seen[curr] = true;
        other_smallest = ((curr < curr_cycle[other_smallest]) && (curr != 0)) ? c : other_smallest;
        curr = PermutationTable[curr];
        c++;
      }
      other_smallest = ((other_smallest == 0) ? int(curr_cycle.size()) : other_smallest);
      for (sint16 i = int(curr_cycle.size()) - 1; i >= 1; i--)
      {
        os << "("
           << ((i >= other_smallest)
                   ? curr_cycle[i] + Rem(CoefficientTable[0] + 1, e) * n
                   : curr_cycle[i] + CoefficientTable[0] * n)
           << ", "
           << ((i >= other_smallest + 1)
                   ? curr_cycle[i - 1] + Rem(CoefficientTable[0] + 1, e) * n
                   : curr_cycle[i - 1] + CoefficientTable[0] * n)
           << ") ";
      }
      os << curr_cycle[0] + CoefficientTable[0] * n << " ";
    }

    for (sint16 i = 1; i <= n; ++i)
    {
      if (!seen[i])
      {
        seen[i] = true;
        c = 0;
        other_smallest = 0;
        cycle_type = CoefficientTable[i];
        if (CoefficientTable[i] == e - 1)
        {
          other_smallest = c + 1;
        }
        curr = PermutationTable[i];
        curr_cycle.clear();
        curr_cycle.push_back(i);
        while (curr != i)
        {
          curr_cycle.push_back(curr);
          seen[curr] = true;
          c++;
          if (CoefficientTable[curr] == e - 1)
          {
            other_smallest = c + 1;
          }
          cycle_type += CoefficientTable[curr];
          curr = PermutationTable[curr];
        }
        if (cycle_type == 2)
        {
          exit(1);
        };
        // If cycle_type == 0, then the cycle is short.
        if (Rem(cycle_type, e) == 0)
        {
          for (sint16 i = int(curr_cycle.size()) - 1; i >= 1; i--)
          {
            os << "("
               << ((other_smallest + i >= int(curr_cycle.size()))
                       ? curr_cycle[other_smallest + i - int(curr_cycle.size())] + n
                       : curr_cycle[other_smallest + i])
               << ", "
               << ((other_smallest + i - 1 >= int(curr_cycle.size()))
                       ? curr_cycle[other_smallest + i - 1 - int(curr_cycle.size())] + n
                       : curr_cycle[other_smallest + i - 1])
               << ") ";
          }
        }
        // Otherwise, it is long.
        else
        {
          other_smallest = ((other_smallest == 0) ? int(curr_cycle.size()) : other_smallest);
          for (sint16 i = int(curr_cycle.size()) - 1; i >= 1; i--)
          {
            os << "("
               << ((i >= other_smallest)
                       ? curr_cycle[i] + n
                       : curr_cycle[i])
               << ", "
               << ((i >= other_smallest + 1)
                       ? curr_cycle[i - 1] + n
                       : curr_cycle[i - 1])
               << ") ";
          }
          os << curr_cycle[0] << " "
             << curr_cycle[0] + n << " ";
        }
      }
    }
  }

  void ComplexDualBraidUnderlying::OfString(std::string str)
  {
    sint16 n = GetParameter().n, e = GetParameter().e;
    sint16 i, j, q, r;
    size_t pos;
    if (str[0] != '(')
    {
      i = Rem(std::stoi(str) - 1, e * n) + 1;
      q = Rem(Quot(i - 1, n), e);
      r = Rem(i - 1, n) + 1;
      Identity();
      PermutationTable[0] = r;
      PermutationTable[r] = 0;
      CoefficientTable[0] = q;
      CoefficientTable[r] = Rem(-q, e);
      return;
    }
    str.erase(0, 1);
    sint16 l = str.length();
    if (str[l - 1] != ')')
    {
      throw InvalidStringError();
    }
    str.erase(l - 1);
    pos = str.find(',');
    i = std::stoi(str.substr(0, pos));
    j = std::stoi(str.substr(pos + 1, str.length()));
    if ((std::abs(i - j) < n) && (i != j))
    {
      sint16 m = std::min(i, j);
      j = std::max(i, j);
      i = m;
      i = Rem(i - 1, n) + 1;
      j = Rem(j - 1, n) + 1;
      j = (j <= i) ? j + n : j;
      Identity();
      PermutationTable[i] = (j > n) ? j - n : j;
      PermutationTable[(j > n) ? j - n : j] = i;
      CoefficientTable[i] = (j > n) ? 1 : 0;
      CoefficientTable[(j > n) ? j - n : j] = (j > n) ? e - 1 : 0;
    }
    else
    {
      throw InvalidStringError();
    }
  }

  void ComplexDualBraidUnderlying::AssignPartition(sint16 *x) const
  {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::vector<sint16> curr_cycle;
    sint16 other_smallest, cycle_type, curr, c;
    for (sint16 i = 1; i <= e * n; ++i)
    {
      x[i] = -1;
    }
    x[0] = 0;

    curr = PermutationTable[0];
    if (curr != 0)
    {
      other_smallest = 0;
      c = 0;
      while (curr != 0)
      {
        curr_cycle.push_back(curr);
        other_smallest = ((curr < curr_cycle[other_smallest]) && (curr != 0)) ? c : other_smallest;
        curr = PermutationTable[curr];
        c++;
      }
      if (other_smallest != 0)
      {
        for (sint16 k = 0; k < other_smallest; k++)
        {
          for (sint16 l = 0; l < e - 1; l++)
          {
            x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
          }
          x[curr_cycle[k] + (e - 1) * n] = curr_cycle[other_smallest];
          x[curr_cycle[k] + Rem(CoefficientTable[0], e) * n] = 0;
        }
        for (sint16 k = other_smallest; k < int(curr_cycle.size()); k++)
        {
          for (sint16 l = 0; l < e - 1; l++)
          {
            x[curr_cycle[k] + (l + 1) * n] = curr_cycle[0] + l * n;
          }
          x[curr_cycle[k]] = curr_cycle[other_smallest];
          x[curr_cycle[k] + Rem(CoefficientTable[0] + 1, e) * n] = 0;
        }
      }
      else
      {
        for (sint16 k = 0; k < int(curr_cycle.size()); k++)
        {
          for (sint16 l = 0; l < e; l++)
          {
            x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
          }
          x[curr_cycle[k] + Rem(CoefficientTable[0], e) * n] = 0;
        }
      }
    }

    for (sint16 i = 1; i <= n; ++i)
    {
      if (x[i] < 0)
      {
        c = 0;
        other_smallest = 0;
        cycle_type = CoefficientTable[i];
        if (CoefficientTable[i] == e - 1)
        {
          other_smallest = 1;
        }
        curr = PermutationTable[i];
        curr_cycle.clear();
        curr_cycle.push_back(i);
        while (curr != i)
        {
          curr_cycle.push_back(curr);
          c++;
          if (CoefficientTable[curr] == e - 1)
          {
            other_smallest = c + 1;
          }
          cycle_type += CoefficientTable[curr];
          curr = PermutationTable[curr];
        }
        // If cycle_type == 0, then the cycle is short.
        if (Rem(cycle_type, e) == 0)
        {
          if (other_smallest != 0)
          {
            for (sint16 k = 0; k < other_smallest; k++)
            {
              x[curr_cycle[k]] = i;
              for (sint16 l = 0; l < e - 1; l++)
              {
                x[curr_cycle[k] + (l + 1) * n] = curr_cycle[other_smallest] + l * n;
              }
            }
            for (sint16 k = other_smallest; k < int(curr_cycle.size()); k++)
            {
              x[curr_cycle[k] + (e - 1) * n] = i;
              for (sint16 l = 0; l < e - 1; l++)
              {
                x[curr_cycle[k] + l * n] = curr_cycle[other_smallest] + l * n;
              }
            }
          }
          else
          {
            for (sint16 k = 0; k < int(curr_cycle.size()); k++)
            {
              for (sint16 l = 0; l < e; l++)
              {
                x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
              }
            }
          }
        }
        // Otherwise, it is long.
        else
        {
          for (sint16 k = 0; k < int(curr_cycle.size()); k++)
          {
            for (sint16 l = 0; l < e; l++)
            {
              x[curr_cycle[k] + l * n] = 0;
            }
          }
        }
      }
    }

    GetParameter().check_non_crossing(x);
  }

  // We assume e > 1 (use the braid group implementation for e = 1).
  void ComplexDualBraidUnderlying::OfPartition(const sint16 *x)
  {
    thread_local sint16 z[MaxBraidIndex + 1];
    sint16 min_cycle_0 = 0, max_cycle_0 = 0;
    sint16 n = GetParameter().n, e = GetParameter().e, r;

    for (sint16 i = 0; i <= n; ++i)
    {
      z[i] = -1;
      PermutationTable[i] = -1;
      CoefficientTable[i] = -1;
    }
    // First find short symmetrical cycles.
    for (sint16 i = 2 * n - 1; i >= n + 1; --i)
    {
      if ((x[i] <= n) && (x[i] >= 1))
      {
        r = i - n;
        if (z[x[i]] == -1)
        {
          PermutationTable[r] = x[i];
          CoefficientTable[r] = e - 1;
          z[x[i]] = r;
        }
        else
        {
          PermutationTable[r] = z[x[i]];
          CoefficientTable[r] = 0;
          z[x[i]] = r;
        }
      }
    }
    for (sint16 i = n; i >= 1; --i)
    {
      if ((x[i] <= n) && (x[i] >= 1) && (x[i + n] > n))
      {
        if ((z[x[i]] == -1))
        {
          PermutationTable[i] = x[i];
          CoefficientTable[i] = 0;
        }
        else
        {
          CoefficientTable[i] = (z[x[i]] < i) ? 1 : 0;
          PermutationTable[i] = z[x[i]];
        }
        z[x[i]] = i;
      }
    }

    // Then consider the part with 0 in it.
    for (sint16 i = 1; i <= e * n; i++)
    {
      if (x[i] == 0)
      {
        min_cycle_0 = i;
        break;
      }
    }
    if (min_cycle_0 != 0)
    {
      // Determine if it is long symmetric.
      if (x[Rem(min_cycle_0 + n - 1, e * n) + 1] == 0)
      {
        CoefficientTable[0] = e - 1;
        PermutationTable[0] = 0;
        for (sint16 i = n; i >= 1; i--)
        {
          if (x[i] == 0)
          {
            if (z[x[i]] == -1)
            {
              PermutationTable[i] = min_cycle_0;
              CoefficientTable[i] = 1;
            }
            else
            {
              PermutationTable[i] = z[x[i]];
              CoefficientTable[i] = 0;
            }
            z[x[i]] = i;
          }
        }
      }
      // Assymetric case.
      else
      {
        for (sint16 i = e * n; i >= 1; i--)
        {
          if (x[i] == 0)
          {
            max_cycle_0 = i;
            break;
          }
        }
        if ((min_cycle_0 <= n) && (max_cycle_0 > (e - 1) * n))
        {
          for (sint16 i = (e - 1) * n + 1; i <= e * n; i++)
          {
            if (x[i] == 0)
            {
              min_cycle_0 = i;
              break;
            }
          }
          for (sint16 i = n; i >= 1; i--)
          {
            if (x[i] == 0)
            {
              max_cycle_0 = i;
              break;
            }
          }
        }
        sint16 q_min = Quot(min_cycle_0 - 1, n), q_max = Quot(max_cycle_0 - 1, n);
        sint16 r_min = Rem(min_cycle_0 - 1, n) + 1;
        z[0] = 0;
        CoefficientTable[0] = q_min;
        PermutationTable[0] = r_min;

        for (sint16 i = n - 1; i >= n - r_min + 1; --i)
        {
          sint16 i_en = Rem(i + min_cycle_0 - 1, e * n) + 1;
          r = Rem(i + min_cycle_0 - 1, n) + 1;
          if (x[i_en] == 0)
          {
            PermutationTable[r] = z[0];
            CoefficientTable[r] = (z[0] == 0) ? Rem(e - q_max, e) : 0;
            z[0] = r;
          }
        }
        for (sint16 i = n - r_min; i >= 0; --i)
        {
          sint16 i_en = Rem(i + min_cycle_0 - 1, e * n) + 1;
          r = Rem(i + min_cycle_0 - 1, n) + 1;
          if ((x[i_en] == 0))
          {
            if (z[0] == 0)
            {
              PermutationTable[r] = z[0];
              CoefficientTable[r] = Rem(e - q_min, e);
            }
            else
            {
              CoefficientTable[r] = (z[0] < r) ? 1 : 0;
              PermutationTable[r] = z[0];
            }
            z[0] = r;
          }
        }
      }
    }

    // Short symmetric case.
    else
    {
      CoefficientTable[0] = 0;
      PermutationTable[0] = 0;
    }
  }

  ComplexDualBraidUnderlying ComplexDualBraidUnderlying::LeftMeet(const ComplexDualBraidUnderlying &b) const
  {
    thread_local sint16 x[MaxE * MaxBraidIndex + 1], y[MaxE * MaxBraidIndex + 1], z[MaxE * MaxBraidIndex + 1];

    AssignPartition(x);
    b.AssignPartition(y);

    thread_local sint16 P[MaxE * MaxBraidIndex + 1][MaxE * MaxBraidIndex + 1];

    for (sint16 i = GetParameter().e * GetParameter().n; i >= 0; i--)
    {
      P[x[i]][y[i]] = i;
    }

    for (sint16 i = 0; i <= GetParameter().e * GetParameter().n; i++)
    {
      z[i] = P[x[i]][y[i]];
    }

    ComplexDualBraidUnderlying c = ComplexDualBraidUnderlying(*this);

    c.OfPartition(z);

    return c;
  }

  void ComplexDualBraidUnderlying::Identity()
  {
    for (sint16 i = 0; i <= GetParameter().n; i++)
    {
      PermutationTable[i] = i;
      CoefficientTable[i] = 0;
    }
  }

  void ComplexDualBraidUnderlying::Delta()
  {
    sint16 i, n = GetParameter().n;
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = i + 1;
      CoefficientTable[i] = 0;
    }
    PermutationTable[0] = 0;
    CoefficientTable[0] = GetParameter().e - 1;
    PermutationTable[n] = 1;
    CoefficientTable[n] = 1;
  }

  bool ComplexDualBraidUnderlying::Compare(const ComplexDualBraidUnderlying &b) const
  {
    sint16 i;
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter().n, b.GetParameter().n);
    }
    for (i = 0; i <= GetParameter().n; i++)
    {
      if ((PermutationTable[i] != b.PermutationTable[i]) || (CoefficientTable[i] != b.CoefficientTable[i]))
      {
        return false;
      }
    }
    return true;
  };

  ComplexDualBraidUnderlying ComplexDualBraidUnderlying::Inverse() const
  {
    ComplexDualBraidUnderlying f = ComplexDualBraidUnderlying(GetParameter());
    Debug(std::cout);
    std::cout << std::endl;
    sint16 i, n = GetParameter().n;
    for (i = 0; i <= n; i++)
    {
      if ((PermutationTable[i] > n) || (PermutationTable[i] < 0))
      {
        std::cout << PermutationTable[i] << std::endl;
        exit(1);
      }
      f.PermutationTable[PermutationTable[i]] = i;
      f.CoefficientTable[PermutationTable[i]] = CoefficientTable[i] == 0 ? 0 : n - CoefficientTable[i] - 1;
    }
    return f;
  };

  ComplexDualBraidUnderlying ComplexDualBraidUnderlying::Product(const ComplexDualBraidUnderlying &b) const
  {
    ComplexDualBraidUnderlying f = ComplexDualBraidUnderlying(GetParameter());
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter().n, b.GetParameter().n);
    }
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 0; i <= n; i++)
    {
      f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
      f.CoefficientTable[i] = Rem(b.CoefficientTable[PermutationTable[i]] + CoefficientTable[i], e);
    }
    return f;
  };

  ComplexDualBraidUnderlying ComplexDualBraidUnderlying::LeftComplement(const ComplexDualBraidUnderlying &b) const
  {
    return b.Product(Inverse());
  };

  ComplexDualBraidUnderlying ComplexDualBraidUnderlying::RightComplement(const ComplexDualBraidUnderlying &b) const
  {
    return Inverse().Product(b);
  };

  void ComplexDualBraidUnderlying::DeltaConjugate(sint16 k)
  {
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    sint16 q = Quot(k, n), r = Rem(k, n), q_e = Rem(q, e);

    ComplexDualBraidUnderlying delta_k = ComplexDualBraidUnderlying(GetParameter());

    delta_k.PermutationTable[0] = 0;
    delta_k.CoefficientTable[0] = Rem(-k, e);
    for (i = 1; i <= n - r; i++)
    {
      delta_k.PermutationTable[i] = Rem(i + k - 1, n) + 1;
      delta_k.CoefficientTable[i] = q_e;
    }
    q_e += 1;
    q_e = q_e == e ? 0 : q_e;
    for (i = n - r + 1; i <= n; i++)
    {
      delta_k.PermutationTable[i] = Rem(i + k - 1, n) + 1;
      delta_k.CoefficientTable[i] = q_e;
    }

    *this = delta_k.Inverse().Product((*this).Product(delta_k));
  }

  void ComplexDualBraidUnderlying::Randomize()
  {
    throw NonRandomizable();
  }

  std::vector<ComplexDualBraidUnderlying> ComplexDualBraidUnderlying::Atoms() const
  {
    ComplexDualBraidParameter p = GetParameter();
    std::vector<ComplexDualBraidUnderlying> atoms;
    ComplexDualBraidUnderlying atom = ComplexDualBraidUnderlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 1; i <= n; i++)
    {
      for (sint16 j = i + 1; j <= n; j++)
      {
        atom.Identity();
        atom.PermutationTable[i] = j;
        atom.PermutationTable[j] = i;
        atoms.push_back(atom);
      }
      for (sint16 j = 1; j < i; j++)
      {
        atom.Identity();
        atom.PermutationTable[i] = j;
        atom.PermutationTable[j] = i;
        atom.CoefficientTable[i] = 1;
        atom.CoefficientTable[j] = -1 + e;
        atoms.push_back(atom);
      }
    }
    for (sint16 k = 0; k < e; k++)
    {
      for (sint16 i = 1; i <= n; i++)
      {
        atom.Identity();
        atom.PermutationTable[0] = i;
        atom.PermutationTable[i] = 0;
        atom.CoefficientTable[0] = k;
        atom.CoefficientTable[i] = -k + e;
        atoms.push_back(atom);
      }
    }
    return atoms;
  }

  std::ostream &operator<<(std::ostream &os, const ComplexDualBraidUnderlying &u)
  {
    u.Print(os);
    return os;
  }
}
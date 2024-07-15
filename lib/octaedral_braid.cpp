#include "octaedral_braid.h"
#include <string>
#include <iostream>

namespace CGarside
{

  sint16 BDualBraidUnderlying::GetParameter() const
  {
    return PresentationParameter;
  }

  sint16 BDualBraidUnderlying::LatticeHeight() const
  {
    return GetParameter();
  }

  BDualBraidUnderlying &BDualBraidUnderlying::Assign(const BDualBraidUnderlying &a)
  {
    if (&a != this)
    {
      for (sint16 i = 1; i <= 2 * GetParameter(); ++i)
      {
        PermutationTable[i] = a.PermutationTable[i];
      }
    }
    return *this;
  };

  BDualBraidUnderlying &BDualBraidUnderlying::operator=(const BDualBraidUnderlying &a)
  {
    return Assign(a);
  }

  BDualBraidUnderlying::BDualBraidUnderlying(sint16 n)
      : PresentationParameter(n)
  {
    PermutationTable = new sint16[2 * n + 1];
  }

  BDualBraidUnderlying::BDualBraidUnderlying(const BDualBraidUnderlying &a)
      : PresentationParameter(a.GetParameter())
  {
    PermutationTable = new sint16[2 * a.GetParameter() + 1];
    sint16 i;
    for (i = 1; i <= 2 * a.GetParameter(); i++)
    {
      PermutationTable[i] = a.PermutationTable[i];
    }
  }

  void BDualBraidUnderlying::Print(std::ostream &os) const
  {
    // Recall that a band braid is represented by decreasing cycles.
    sint16 i, j, n = GetParameter();
    std::vector<sint16> curr_cycle;
    std::vector<bool> seen(n + 1, false);
    for (i = 1; i <= n; i++)
    {
      seen[i] = false;
    }
    for (i = 1; i <= n; ++i)
    {
      if (!seen[i])
      {
        curr_cycle.clear();
        j = i;
        while (j < PermutationTable[j])
        {
          curr_cycle.push_back(j);
          seen[j] = true;
          seen[((j + n - 1) % (2 * n)) + 1] = true;
          j = PermutationTable[j];
        }
        curr_cycle.push_back(j);
        seen[j] = true;
        seen[((j + n - 1) % (2 * n)) + 1] = true;
        for (sint16 l = int(curr_cycle.size()); l >= 1; --l)
        {
          if (l == ((j + n - 1) % (2 * n)) + 1)
          {
            os << j << " ";
            break;
          }
          os << "(" << curr_cycle[l] << ", " << curr_cycle[l - 1] << ") ";
        }
      }
    }
  }

  void BDualBraidUnderlying::OfString(std::string &str)
  {
    sint16 n = GetParameter();
    sint16 i, j, k;
    size_t pos;
    if (str[0] != '(')
    {
      i = std::stoi(str);
      if (i <= n)
      {
        for (k = 1; k <= 2 * n; k++)
        {
          PermutationTable[k] = k;
        }
        j = ((i + n - 1) % (2 * n)) + 1;
        PermutationTable[i] = j;
        PermutationTable[j] = i;
        return;
      }
      else
      {
        throw InvalidStringError();
      }
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
    if (i <= 2 * n && j <= 2 * n)
    {
      for (k = 1; k <= 2 * n; k++)
      {
        PermutationTable[k] = k;
      }
      PermutationTable[i] = j;
      PermutationTable[j] = i;
      PermutationTable[((i + n - 1) % (2 * n)) + 1] = ((j + n - 1) % (2 * n)) + 1;
      PermutationTable[((j + n - 1) % (2 * n)) + 1] = ((i + n - 1) % (2 * n)) + 1;
    }
    else
    {
      throw InvalidStringError();
    }
  }

  void BDualBraidUnderlying::AssignDCDT(sint16 *x) const
  {
    for (sint16 i = 1; i <= 2 * GetParameter(); ++i)
      x[i] = 0;
    for (sint16 i = 1; i <= 2 * GetParameter(); ++i)
    {
      if (x[i] == 0)
        x[i] = i;
      if (PermutationTable[i] > i)
        x[PermutationTable[i]] = x[i];
    }
  }

  void BDualBraidUnderlying::OfDCDT(const sint16 *x)
  {
    thread_local sint16 z[2 * MaxBraidIndex + 1];

    for (sint16 i = 1; i <= 2 * GetParameter(); ++i)
      z[i] = 0;
    for (sint16 i = 2 * GetParameter(); i >= 1; --i)
    {
      PermutationTable[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
      z[x[i]] = i;
    }
  }

  BDualBraidUnderlying BDualBraidUnderlying::LeftMeet(const BDualBraidUnderlying &b) const
  {
    thread_local sint16 x[2 * MaxBraidIndex + 1], y[2 * MaxBraidIndex + 1], z[2 * MaxBraidIndex + 1];

    AssignDCDT(x);
    b.AssignDCDT(y);

    thread_local sint16 P[2 * MaxBraidIndex + 1][2 * MaxBraidIndex + 1];

    for (sint16 i = 2 * GetParameter(); i >= 1; i--)
    {
      P[x[i]][y[i]] = i;
    }

    for (sint16 i = 1; i <= 2 * GetParameter(); i++)
    {
      z[i] = P[x[i]][y[i]];
    }

    BDualBraidUnderlying c = BDualBraidUnderlying(*this);

    c.OfDCDT(z);

    return c;
  }

  BDualBraidUnderlying BDualBraidUnderlying::RightMeet(const BDualBraidUnderlying &b) const
  {
    return LeftMeet(b);
  }

  void BDualBraidUnderlying::Identity()
  {
    for (sint16 i = 1; i <= 2 * GetParameter(); i++)
    {
      PermutationTable[i] = i;
    }
  }

  void BDualBraidUnderlying::Delta()
  {
    sint16 i, n = GetParameter();
    for (i = 1; i < 2 * n; i++)
    {
      PermutationTable[i] = i + 1;
    }
    PermutationTable[2 * n] = 1;
  }

  bool BDualBraidUnderlying::Compare(const BDualBraidUnderlying &b) const
  {
    sint16 i;
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter(), b.GetParameter());
    }
    for (i = 1; i <= GetParameter(); i++)
    {
      if (PermutationTable[i] != b.PermutationTable[i])
      {
        return false;
      }
    }
    return true;
  };

  BDualBraidUnderlying BDualBraidUnderlying::Inverse() const
  {
    BDualBraidUnderlying f = BDualBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= 2 * GetParameter(); i++)
    {
      f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
  };

  BDualBraidUnderlying BDualBraidUnderlying::Product(const BDualBraidUnderlying &b) const
  {
    BDualBraidUnderlying f = BDualBraidUnderlying(GetParameter());
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter(), b.GetParameter());
    }
    sint16 i;
    for (i = 1; i <= 2 * GetParameter(); i++)
    {
      f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
    }
    return f;
  };

  BDualBraidUnderlying BDualBraidUnderlying::LeftComplement(const BDualBraidUnderlying &b) const
  {
    return b.Product(Inverse());
  };

  BDualBraidUnderlying BDualBraidUnderlying::RightComplement(const BDualBraidUnderlying &b) const
  {
    return Inverse().Product(b);
  };

  BDualBraidUnderlying BDualBraidUnderlying::DeltaConjugate(sint16 k) const
  {
    BDualBraidUnderlying under = BDualBraidUnderlying(*this);
    sint16 i, n = GetParameter();

    if (k < 0)
    {
      k = k - 2 * n * k;
    }

    for (i = 1; i < 2 * n; i++)
    {
      under.PermutationTable[i] = (PermutationTable[(i + k - 1) % (2 * n) + 1] - k - 1) % (2 * n) + 1;
    }
    return under;
  }

  void BDualBraidUnderlying::Randomize()
  {
    throw NonRandomizable();
  }

  std::vector<BDualBraidUnderlying> BDualBraidUnderlying::Atoms() const
  {
    sint16 n = GetParameter();
    std::vector<BDualBraidUnderlying> atoms;
    for (sint16 i = 1; i <= n; i++)
    {
      for (sint16 j = i + 1; j <= n; j++)
      {
        BDualBraidUnderlying atom = BDualBraidUnderlying(n);
        atom.Identity();
        atom.PermutationTable[i] = j;
        atom.PermutationTable[j] = i;
        atom.PermutationTable[i + n] = j + n;
        atom.PermutationTable[j + n] = i + n;
        atoms.push_back(atom);
      }
      for (sint16 j = n + i + 1; j <= 2 * n; j++)
      {
        BDualBraidUnderlying atom = BDualBraidUnderlying(n);
        atom.Identity();
        atom.PermutationTable[i] = j;
        atom.PermutationTable[j] = i;
        atom.PermutationTable[i + n] = j - n;
        atom.PermutationTable[j - n] = i + n;
        atoms.push_back(atom);
      }
        BDualBraidUnderlying atom = BDualBraidUnderlying(n);
        atom.Identity();
        atom.PermutationTable[i] = n + i;
        atom.PermutationTable[n + i] = i;
        atoms.push_back(atom);
    }
    return atoms;
  }

}
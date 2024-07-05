#include "band_braid.h"
#include "cbraid.h"
#include <string>
#include <iostream>

namespace CGarside
{

  sint16 BandBraidUnderlying::GetParameter() const
  {
    return PresentationParameter;
  }

  sint16 BandBraidUnderlying::LatticeHeight() const
  {
    return GetParameter();
  }

  BandBraidUnderlying &BandBraidUnderlying::Assign(const BandBraidUnderlying &a)
  {
    if (&a != this)
    {
      for (sint16 i = 1; i <= GetParameter(); ++i)
      {
        PermutationTable[i] = a.PermutationTable[i];
      }
    }
    return *this;
  };

  BandBraidUnderlying &BandBraidUnderlying::operator=(const BandBraidUnderlying &a)
  {
    return Assign(a);
  }

  BandBraidUnderlying::BandBraidUnderlying(sint16 n)
      : PresentationParameter(n)
  {
    PermutationTable = new sint16[n + 1];
  }

  BandBraidUnderlying::BandBraidUnderlying(const BandBraidUnderlying &a)
      : PresentationParameter(a.GetParameter())
  {
    PermutationTable = new sint16[a.GetParameter() + 1];
    sint16 i;
    for (i = 1; i <= a.GetParameter(); i++)
    {
      PermutationTable[i] = a.PermutationTable[i];
    }
  }

  void BandBraidUnderlying::Print(std::ostream &os) const
  {
    // Recall that a band braid is represented by decreasing cycles.
    sint16 i, j, k, n = GetParameter();
    sint16 curr_cycle[n];
    bool seen[n + 1];
    for (i = 1; i <= n; i++)
    {
      seen[i] = false;
    }
    for (i = 1; i <= n; ++i)
    {
      if (!seen[i])
      {
        k = 0;
        j = i;
        while (j < PermutationTable[j])
        {
          curr_cycle[k] = j;
          k++;
          seen[j] = true;
          j = PermutationTable[j];
        }
        curr_cycle[k] = j;
        seen[j] = true;
        for (j = k; j >= 1; --j)
        {
          os << "(" << curr_cycle[j] << ", " << curr_cycle[j - 1] << ")" << " ";
        }
      }
    }
  }

  void BandBraidUnderlying::OfString(std::string &str)
  {
    sint16 n = GetParameter();
    sint16 i, j, k;
    size_t pos;
    if (str[0] != '(')
    {
      throw InvalidStringError();
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
    if (i <= n && j < i)
    {
      for (k = 1; k <= n; k++)
      {
        PermutationTable[k] = k;
      }
      PermutationTable[i] = j;
      PermutationTable[j] = i;
    }
    else
    {
      throw InvalidStringError();
    }
  }

  void BandBraidUnderlying::AssignDCDT(sint16 *x) const
  {
    for (sint16 i = 1; i <= GetParameter(); ++i)
      x[i] = 0;
    for (sint16 i = 1; i <= GetParameter(); ++i)
    {
      if (x[i] == 0)
        x[i] = i;
      if (PermutationTable[i] > i)
        x[PermutationTable[i]] = x[i];
    }
  }

  void BandBraidUnderlying::OfDCDT(const sint16 *x)
  {
    thread_local sint16 z[MaxBraidIndex];

    for (sint16 i = 1; i <= GetParameter(); ++i)
      z[i] = 0;
    for (sint16 i = GetParameter(); i >= 1; --i)
    {
      PermutationTable[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
      z[x[i]] = i;
    }
  }

  BandBraidUnderlying BandBraidUnderlying::LeftMeet(const BandBraidUnderlying &b) const
  {
    thread_local sint16 x[MaxBraidIndex], y[MaxBraidIndex], z[MaxBraidIndex];

    AssignDCDT(x);
    b.AssignDCDT(y);

    thread_local sint16 P[MaxBraidIndex][MaxBraidIndex];

    for (sint16 i = GetParameter(); i >= 1; i--)
    {
      P[x[i]][y[i]] = i;
    }

    for (sint16 i = 1; i <= GetParameter(); i++)
    {
      z[i] = P[x[i]][y[i]];
    }

    BandBraidUnderlying c = BandBraidUnderlying(*this);

    c.OfDCDT(z);

    return c;
  }

  BandBraidUnderlying BandBraidUnderlying::RightMeet(const BandBraidUnderlying &b) const
  {
    return LeftMeet(b);
  }

  void BandBraidUnderlying::Identity()
  {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = i;
    }
  }

  void BandBraidUnderlying::Delta()
  {
    sint16 i, n = GetParameter();
    for (i = 1; i < n; i++)
    {
      PermutationTable[i] = i + 1;
    }
    PermutationTable[n] = 1;
  }

  bool BandBraidUnderlying::Compare(const BandBraidUnderlying &b) const
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

  BandBraidUnderlying BandBraidUnderlying::Inverse() const
  {
    BandBraidUnderlying f = BandBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++)
    {
      f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
  };

  BandBraidUnderlying BandBraidUnderlying::Product(const BandBraidUnderlying &b) const
  {
    BandBraidUnderlying f = BandBraidUnderlying(GetParameter());
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes(GetParameter(), b.GetParameter());
    }
    sint16 i;
    for (i = 1; i <= GetParameter(); i++)
    {
      f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
    }
    return f;
  };

  BandBraidUnderlying BandBraidUnderlying::LeftComplement(const BandBraidUnderlying &b) const
  {
    return b.Product(Inverse());
  };

  BandBraidUnderlying BandBraidUnderlying::RightComplement(const BandBraidUnderlying &b) const
  {
    return Inverse().Product(b);
  };

  BandBraidUnderlying BandBraidUnderlying::DeltaConjugate(sint16 k) const
  {
    BandBraidUnderlying under = BandBraidUnderlying(*this);
    sint16 i, n = GetParameter();

    if (k < 0)
    {
      k = k - n * k;
    }

    for (i = 1; i < n; i++)
    {
      under.PermutationTable[i] = (PermutationTable[(i + k - 1) % n + 1] - k - 1) % n + 1;
    }
    return under;
  }

  void BandBraidUnderlying::Randomize()
  {
    CBraid::BandPresentation pres = CBraid::BandPresentation(GetParameter());
    pres.Randomize(PermutationTable);
  }

  std::vector<BandBraidUnderlying> BandBraidUnderlying::Atoms() const
  {
    sint16 n = GetParameter();
    std::vector<BandBraidUnderlying> atoms;
    for (sint16 i = 1; i <= n; i++)
    {
      for (sint16 j = 1; j < i; j++)
      {
        BandBraidUnderlying atom = BandBraidUnderlying(n);
        atom.Identity();
        atom.PermutationTable[i] = j;
        atom.PermutationTable[j] = i;
        atoms.push_back(atom);
      }
    }
    return atoms;
  }

}

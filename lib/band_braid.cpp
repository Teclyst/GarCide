#include "band_braid.h"
#include "cbraid.h"
#include <string>
#include <iostream>

namespace CGarside
{

  sint16 BandBraidUnderlying::Index() const
  {
    return PresentationIndex;
  }

  BandBraidUnderlying &BandBraidUnderlying::Assign(const BandBraidUnderlying &a)
  {
    if (&a != this)
    {
      for (sint16 i = 1; i <= Index(); ++i)
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
      : PresentationIndex(n)
  {
    PermutationTable = new sint16[n + 1];
  }

  BandBraidUnderlying::BandBraidUnderlying(const BandBraidUnderlying &a)
      : PresentationIndex(a.Index())
  {
    PermutationTable = new sint16[a.Index() + 1];
    sint16 i;
    for (i = 1; i <= a.Index(); i++)
    {
      PermutationTable[i] = a.PermutationTable[i];
    }
  }

  void BandBraidUnderlying::Print(std::ostream &os)
  {
    // Recall that a band braid is represented by decreasing cycles.
    sint16 i, j, k, n = Index();
    sint16 curr_cycle[n];
    bool seen[n + 1];
    for (i = 1; i <= n; i++)
    {
      seen[i] = false;
    }
    for (i = 1; i <= n; ++i)
    {
      if (not seen[i])
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
    sint16 n = Index();
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
    for (sint16 i = 1; i <= Index(); ++i)
      x[i] = 0;
    for (sint16 i = 1; i <= Index(); ++i)
    {
      if (x[i] == 0)
        x[i] = i;
      if (PermutationTable[i] > i)
        x[PermutationTable[i]] = x[i];
    }
  }

  void BandBraidUnderlying::OfDCDT(const sint16 *x)
  {
    static sint16 z[MaxBraidIndex];

    for (sint16 i = 1; i <= Index(); ++i)
      z[i] = 0;
    for (sint16 i = Index(); i >= 1; --i)
    {
      PermutationTable[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
      z[x[i]] = i;
    }
  }

  BandBraidUnderlying BandBraidUnderlying::LeftMeet(const BandBraidUnderlying &b) const
  {
    static sint16 x[MaxBraidIndex], y[MaxBraidIndex], z[MaxBraidIndex];

    AssignDCDT(x);
    b.AssignDCDT(y);

    static sint16 P[MaxBraidIndex][MaxBraidIndex];

    for (sint16 i = Index(); i >= 1; i--)
    {
      P[x[i]][y[i]] = i;
    }

    for (sint16 i = 1; i <= Index(); i++)
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
    sint16 i, n = Index();
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = i;
    }
  }

  void BandBraidUnderlying::Delta()
  {
    sint16 i, n = Index();
    for (i = 1; i < n; i++)
    {
      PermutationTable[i] = i + 1;
    }
    PermutationTable[n] = 1;
  }

  bool BandBraidUnderlying::Compare(const BandBraidUnderlying &b) const
  {
    sint16 i;
    if (Index() != b.Index())
    {
      throw NonMatchingIndexes();
    }
    for (i = 1; i <= Index(); i++)
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
    BandBraidUnderlying f = BandBraidUnderlying(Index());
    sint16 i;
    for (i = 1; i <= Index(); i++)
    {
      f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
  };

  BandBraidUnderlying BandBraidUnderlying::Product(const BandBraidUnderlying &b) const
  {
    BandBraidUnderlying f = BandBraidUnderlying(Index());
    if (Index() != b.Index())
    {
      throw NonMatchingIndexes();
    }
    sint16 i;
    for (i = 1; i <= Index(); i++)
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
    sint16 i, n = Index();

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

  void BandBraidUnderlying::Randomize() {
    CBraid::BandPresentation pres = CBraid::BandPresentation(Index());
    pres.Randomize(PermutationTable);
  }

}

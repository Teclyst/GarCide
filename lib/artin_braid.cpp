#include "artin_braid.h"

namespace CGarside
{

  void ArtinBraidUnderlying::MeetSub(const sint16 *a, const sint16 *b, sint16 *r,
                                     sint16 s, sint16 t)
  {
    static sint16 u[MaxBraidIndex], v[MaxBraidIndex], w[MaxBraidIndex];

    if (s >= t)
      return;
    sint16 m = (s + t) / 2;
    MeetSub(a, b, r, s, m);
    MeetSub(a, b, r, m + 1, t);

    u[m] = a[r[m]];
    v[m] = b[r[m]];
    if (s < m)
    {
      for (sint16 i = m - 1; i >= s; --i)
      {
        u[i] = std::min(a[r[i]], u[i + 1]);
        v[i] = std::min(b[r[i]], v[i + 1]);
      }
    }
    u[m + 1] = a[r[m + 1]];
    v[m + 1] = b[r[m + 1]];
    if (t > m + 1)
    {
      for (sint16 i = m + 2; i <= t; ++i)
      {
        u[i] = std::max(a[r[i]], u[i - 1]);
        v[i] = std::max(b[r[i]], v[i - 1]);
      }
    }

    sint16 p = s;
    sint16 q = m + 1;
    for (sint16 i = s; i <= t; ++i)
      w[i] = ((p > m) || (q <= t && u[p] > u[q] && v[p] > v[q])) ? r[q++] : r[p++];
    for (sint16 i = s; i <= t; ++i)
      r[i] = w[i];
  }

  sint16 ArtinBraidUnderlying::GetParameter() const
  {
    return PresentationParameter;
  };

  ArtinBraidUnderlying &ArtinBraidUnderlying::Assign(const ArtinBraidUnderlying &a)
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

  ArtinBraidUnderlying &ArtinBraidUnderlying::operator=(const ArtinBraidUnderlying &a)
  {
    return Assign(a);
  }

  ArtinBraidUnderlying::ArtinBraidUnderlying(sint16 n)
      : PresentationParameter(n)
  {
    PermutationTable = new sint16[n + 1];
  }

  ArtinBraidUnderlying::ArtinBraidUnderlying(const ArtinBraidUnderlying &a)
      : PresentationParameter(a.GetParameter())
  {
    PermutationTable = new sint16[a.GetParameter() + 1];
    sint16 i;
    for (i = 1; i <= a.GetParameter(); i++)
    {
      PermutationTable[i] = a.PermutationTable[i];
    }
  }

  void ArtinBraidUnderlying::OfString(const std::string &str)
  {
    sint16 n = GetParameter();
    sint16 i;
    sint16 j;
    try
    {
      i = std::stoi(str);
    }
    catch (std::invalid_argument)
    {
      throw InvalidStringError();
    }
    if (i == n)
    {
      for (j = 1; j <= n; j++)
      {
        PermutationTable[j] = n + 1 - j;
      }
    }
    else if (i < n)
    {
      for (j = 1; j <= n; j++)
      {
        PermutationTable[j] = j;
      }
      PermutationTable[i] = i + 1;
      PermutationTable[i + 1] = i;
    }
    else
    {
      throw InvalidStringError();
    }
  };

  void ArtinBraidUnderlying::Print(std::ostream &os) const
  {
    sint16 i, j, k, n = GetParameter();

    ArtinBraidUnderlying c = ArtinBraidUnderlying(*this); 

    for (i = 2; i <= n; i++)
    {
      for (j = i; j > 1 && c.PermutationTable[j] < c.PermutationTable[j - 1]; j--)
      {
        os << j - 1 << " ";
        k = c.PermutationTable[j];
        c.PermutationTable[j] = c.PermutationTable[j - 1];
        c.PermutationTable[j - 1] = k;
      }
    }
  };

  sint16 ArtinBraidUnderlying::LatticeHeight() const {
    sint16 n = GetParameter();
    return n * (n - 1) / 2;
  }

  void ArtinBraidUnderlying::Identity()
  {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = i;
    }
  };

  void ArtinBraidUnderlying::Delta()
  {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = n + 1 - i;
    }
  };

  ArtinBraidUnderlying ArtinBraidUnderlying::LeftMeet(const ArtinBraidUnderlying &b) const
  {
    static sint16 s[MaxBraidIndex];

    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes();
    }

    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());

    for (sint16 i = 1; i <= GetParameter(); ++i)
      s[i] = i;
    MeetSub(PermutationTable, b.PermutationTable, s, 1, GetParameter());
    for (sint16 i = 1; i <= GetParameter(); ++i)
      f.PermutationTable[s[i]] = i;

    return f;
  };

  ArtinBraidUnderlying ArtinBraidUnderlying::RightMeet(const ArtinBraidUnderlying &b) const
  {
    static sint16 u[MaxBraidIndex], v[MaxBraidIndex];

    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes();
    }

    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());

    for (sint16 i = 1; i <= GetParameter(); ++i)
    {
      u[PermutationTable[i]] = i;
      v[b.PermutationTable[i]] = i;
    }
    for (sint16 i = 1; i <= GetParameter(); ++i)
      f.PermutationTable[i] = i;
    MeetSub(u, v, f.PermutationTable, 1, GetParameter());

    return f;
  };

  bool ArtinBraidUnderlying::Compare(const ArtinBraidUnderlying &b) const
  {
    sint16 i;
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes();
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

  ArtinBraidUnderlying ArtinBraidUnderlying::Inverse() const
  {
    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++)
    {
      f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
  };

  ArtinBraidUnderlying ArtinBraidUnderlying::Product(const ArtinBraidUnderlying &b) const
  {
    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());
    if (GetParameter() != b.GetParameter())
    {
      throw NonMatchingIndexes();
    }
    sint16 i;
    for (i = 1; i <= GetParameter(); i++)
    {
      f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
    }
    return f;
  };

  ArtinBraidUnderlying ArtinBraidUnderlying::LeftComplement(const ArtinBraidUnderlying &b) const
  {
    return b.Product(Inverse());
  };

  ArtinBraidUnderlying ArtinBraidUnderlying::RightComplement(const ArtinBraidUnderlying &b) const
  {
    return Inverse().Product(b);
  };

  void ArtinBraidUnderlying::Randomize()
  {
    for (sint16 i = 1; i <= GetParameter(); ++i)
      PermutationTable[i] = i;
    for (sint16 i = 1; i < GetParameter(); ++i)
    {
      sint16 j = i + sint16(std::rand() / (RAND_MAX + 1.0) * (GetParameter() - i + 1));
      sint16 z = PermutationTable[i];
      PermutationTable[i] = PermutationTable[j];
      PermutationTable[j] = z;
    }
  };

  std::list<ArtinBraidUnderlying> ArtinBraidUnderlying::Atoms() const
  {
    sint16 i, n = GetParameter();
    std::list<ArtinBraidUnderlying> atoms;
    for (i = 1; i <= n - 1; i++)
    {
      ArtinBraidUnderlying atom = ArtinBraidUnderlying(n);
      atom.OfString(std::to_string(i));
      atoms.push_front(atom);
    }
    return atoms;
  };

  ArtinBraidUnderlying ArtinBraidUnderlying::DeltaConjugate(sint16 k) const
  {
    sint16 i, n = GetParameter();
    if (k % 2 == 0)
    {
      return ArtinBraidUnderlying(*this);
    }
    else
    {
      ArtinBraidUnderlying under = ArtinBraidUnderlying(*this);
      for (i = 1; i <= n / 2; i++)
      {
        sint16 u = under.PermutationTable[i];
        under.PermutationTable[i] = under.PermutationTable[n - i + 1];
        under.PermutationTable[n - i + 1] = u;
      }
      return under;
    }
  }

}
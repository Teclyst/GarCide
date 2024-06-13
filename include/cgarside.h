#include <list>
#include <string>
#include <ostream>
#include <stdexcept>

namespace CGarside {

typedef char sint8;
typedef unsigned char uint8;
/* JLT: Don't use short ints... */
// typedef short sint16;
// typedef unsigned short uint16;
typedef int sint16;
typedef unsigned int uint16;
typedef int sint32;
typedef unsigned int uint32;
typedef long long sint64;
typedef unsigned long long uint64;

// Maximum braid index.
const sint16 MaxBraidIndex = 300;

// Exceptions.
struct InvalidStringError
{
};
struct NonMatchingIndexes
{
};

// A class for the underlying objects for canonical factors
// in the Artin presentation braid group case.
// In this case, permutations.
class ArtinBraidUnderlying
{

protected:
  sint16 PresentationIndex;

  sint16 *PermutationTable;

public:
  sint16 Index() const
  {
    return PresentationIndex;
  };

  // Constructor
  ArtinBraidUnderlying(sint16 n)
  {
    PresentationIndex = n;
    PermutationTable = new sint16[n + 1];
  }

  ArtinBraidUnderlying Copy()
  {
    ArtinBraidUnderlying copy = ArtinBraidUnderlying(Index());
    sint16 i;
    for (i = 1; i <= Index(); i++)
    {
      copy.PermutationTable[i] = PermutationTable[i];
    }
    return copy;
  }

  void OfString(std::string str) const
  {
    sint16 n = Index();
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

  // Print to os. Be wary, as it side-effects!
  void Print(std::ostream &os)
  {
    sint16 i, j, k, n = Index();

    for (i = 2; i <= n; i++)
    {
      for (j = i; j > 1 && PermutationTable[j] < PermutationTable[j - 1]; j--)
      {
        os << j - 1 << " ";
        k = PermutationTable[j];
        PermutationTable[j] = PermutationTable[j - 1];
        PermutationTable[j - 1] = k;
      }
    }
  };

  // Returns the neutral element (here the identity).
  void Neutral()
  {
    sint16 i, n = Index();
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = i;
    }
  };

  // Returns delta.
  void Delta()
  {
    sint16 i, n = Index();
    for (i = 1; i <= n; i++)
    {
      PermutationTable[i] = n + 1 - i;
    }
  };

  // Compute the meet r of the two factors a and b.  A factor is
  // given as the associated permutation, which is viewed as a
  // bijection on the set {1,...n} and represented as an array whose
  // i-th entry is the image of i under the inverse of the
  // permutation (this convention is different from that in the
  // AsiaCrypt 2001 paper of the author).  The range of indices is
  // [1,n], not [0,n[.  We use a C style array of size (n+1) to
  // represent an n-permutation (the first entry is not used).

  // We define the left meet of two factors a and b to be the
  // longest factor r such that a=ra' and b=rb' for some factors a'
  // and b'.  This coincides with the convention of the paper of
  // Birman, Ko, and Lee, but different from that of the article of
  // Thurston (in Epstein's book).  Indeed, Thurston's is the
  // "right" meet in our sense.
  ArtinBraidUnderlying LeftMeet(const ArtinBraidUnderlying b) const
  {
    static sint16 s[MaxBraidIndex];

    if (Index() != b.Index())
    {
      throw NonMatchingIndexes();
    }

    ArtinBraidUnderlying f = ArtinBraidUnderlying(Index());

    for (sint16 i = 1; i <= Index(); ++i)
      s[i] = i;
    MeetSub(PermutationTable, b.PermutationTable, s, 1, Index());
    for (sint16 i = 1; i <= Index(); ++i)
      f.PermutationTable[s[i]] = i;
  };

  ArtinBraidUnderlying RightMeet(const ArtinBraidUnderlying b) const
  {
    static sint16 u[MaxBraidIndex], v[MaxBraidIndex];

    if (Index() != b.Index())
    {
      throw NonMatchingIndexes();
    }

    ArtinBraidUnderlying f = ArtinBraidUnderlying(Index());

    for (sint16 i = 1; i <= Index(); ++i)
    {
      u[PermutationTable[i]] = i;
      v[b.PermutationTable[i]] = i;
    }
    for (sint16 i = 1; i <= Index(); ++i)
      f.PermutationTable[i] = i;
    MeetSub(u, v, f.PermutationTable, 1, Index());
  };

  // Equality check.
  // We check wether the underlying permutation table are (pointwise) equal.
  bool Compare(const ArtinBraidUnderlying b) const
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

  // Computes the factor corresponding to the inverse permutation.
  // Used to simplify complement operation.
  ArtinBraidUnderlying Inverse() const
  {
    ArtinBraidUnderlying f = ArtinBraidUnderlying(Index());
    sint16 i;
    for (i = 1; i <= Index(); i++)
    {
      f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
  };

  // Product under the hypothesis that it is still simple.
  ArtinBraidUnderlying Product(const ArtinBraidUnderlying b) const
  {
    ArtinBraidUnderlying f = ArtinBraidUnderlying(Index());
    if (Index() != b.Index())
    {
      throw NonMatchingIndexes();
    }
    sint16 i;
    for (i = 1; i <= Index(); i++)
    {
      f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
    }
  };

  // Under the assumption a <= b, a.LeftComplement(b) computes
  // The factor c such that ac = b.
  ArtinBraidUnderlying LeftComplement(const ArtinBraidUnderlying b) const
  {
    return Inverse().Product(b);
  };

  ArtinBraidUnderlying RightComplement(const ArtinBraidUnderlying b) const
  {
    return b.Product(Inverse());
  };

  // Generate a random factor.
  void Randomize()
  {
    for (sint16 i = 1; i <= Index(); ++i)
      PermutationTable[i] = i;
    for (sint16 i = 1; i < Index(); ++i)
    {
      sint16 j = i + sint16(std::rand() / (RAND_MAX + 1.0) * (Index() - i + 1));
      sint16 z = PermutationTable[i];
      PermutationTable[i] = PermutationTable[j];
      PermutationTable[j] = z;
    }
  };

  // List of atoms.
  std::list<ArtinBraidUnderlying> Atoms() const
  {
    sint16 i, n = Index();
    std::list<ArtinBraidUnderlying> atoms;
    for (i = 1; i <= n - 1; i++)
    {
      ArtinBraidUnderlying atom = ArtinBraidUnderlying(n);
      atom.OfString(std::to_string(i));
      atoms.push_front(atom);
    }
    return atoms;
  };

private:
  // Subroutine called by LeftMeet() and RightMeet()
  static void MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s, sint16 t);
};

template <class U>
class Factor
{

private:
  U Underlying;

public:
  // Constructor
  Factor(U under): Underlying(under)
  {
  };

  Factor Copy() const
  {
    return Factor(Underlying.Copy());
  }

  void OfString(std::string str)
  {
    Underlying.OfString(str);
  };

  void Print(std::ostream &os)
  {
    Underlying.Print(os);
  };

  void Neutral()
  {
    Underlying.Neutral();
  }

  void Delta()
  {
    Underlying.Delta();
  }

  bool Compare(const Factor b) const
  {
    return Underlying.Compare(b.Underlying);
  }

  bool operator==(const Factor b) const
  {
    return Compare(b);
  }

  bool operator!=(const Factor b) const
  {
    return !Compare(b);
  }

  bool IsNeutral() const
  {
    Factor e = Copy();
    e.Neutral();
    return Underlying.Compare(e);
  }

  bool IsDelta() const
  {
    Factor delta = Copy();
    delta.Delta();
    return Underlying.Compare(delta);
  }

  Factor LeftComplement(const Factor b) const
  {
    return Factor(Underlying.LeftComplement(b.Underlying));
  }

  Factor operator/(const Factor b) const
  {
    return b.LeftComplement(this);
  }

  Factor LeftComplement() const
  {
    Factor delta = Copy();
    delta.Delta();
    return LeftComplement(delta);
  }

  Factor operator~() const
  {
    return LeftComplement();
  }

  Factor RightComplement(const Factor b) const
  {
    return Factor(Underlying.RightComplement(b.Underlying));
  }

  Factor RightComplement() const
  {
    Factor delta = Copy();
    delta.Delta();
    return RightComplement(delta);
  }

  // Conjugation by delta, AKA double complementation.
  Factor DeltaConjugate() const
  {
    ~Left(LeftComplement());
  }

  Factor DeltaConjugate(const sint16 k) const
  {
    Factor conjugate = Copy();
    sint16 i;
    for (i = 0; i < k; i++)
    {
      conjugate = conjugate.DeltaConjugate();
    }
    return conjugate;
  }

  Factor LeftMeet(const Factor b) const
  {
    return Factor(Underlying.LeftMeet(b.Underlying));
  }

  Factor operator^(const Factor b) const
  {
    return LeftMeet(b);
  }

  Factor RightMeet(const Factor b) const
  {
    return Factor(Underlying.RightMeet(b.Underlying));
  }

  // a.IsLeftWeighted(b) returns true if a | b is left weighted, or false otherwise.
  bool IsLeftWeighted(const Factor b) const
  {
    RightComplement().LeftMeet(b).IsTrivial();
  }

  // a.IsRightWeighted(b) returns true if a | b is right weighted, or false otherwise.
  bool IsRightWeighted(const Factor b) const
  {
    LeftMeet(b.LeftComplement()).IsTrivial();
  }

  Factor Product(const Factor b) const
  {
    return Factor(Underlying.Product(b.Underlying));
  }

  Factor operator*(const Factor b) const
  {
    return Product(b);
  }

  void Randomize()
  {
    Underlying.Randomize();
  }

  std::list<Factor> Atoms() const
  {
    std::list<U> atoms = Underlying.Atoms();
    typename std::list<U>::iterator atoms_it;
    std::list<Factor> factor_atoms;
    for (auto const &atoms_it : atoms)
    {
      Factor f = Copy();
      f.Set(atoms_it);
      factor_atoms.push_back(f);
    }
    return factor_atoms;
  }
};

template<class U>
void MakeLeftWeighted(Factor<U> &u, Factor<U> &v){
  Factor<U> t = ~u ^ v;
  v = v / t;
  u = u * t;
}

typedef Factor<ArtinBraidUnderlying> ArtinBraidFactor;

}
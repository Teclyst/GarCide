#include <list>
#include <string>
#include <ostream>
#include <stdexcept>

namespace CGarside
{

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

    ArtinBraidUnderlying(ArtinBraidUnderlying &a)
    {
      PresentationIndex = a.Index();
      sint16 i;
      for (i = 1; i <= a.Index(); i++)
      {
        PermutationTable[i] = a.PermutationTable[i];
      }
    }

    void OfString(const std::string &str) const
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

    // Set to the neutral element (here the identity).
    void Neutral()
    {
      sint16 i, n = Index();
      for (i = 1; i <= n; i++)
      {
        PermutationTable[i] = i;
      }
    };

    // Set to delta.
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
    ArtinBraidUnderlying LeftMeet(const ArtinBraidUnderlying &b) const
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

    ArtinBraidUnderlying RightMeet(const ArtinBraidUnderlying &b) const
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
    bool Compare(const ArtinBraidUnderlying &b) const
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
    ArtinBraidUnderlying Product(const ArtinBraidUnderlying &b) const
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
    ArtinBraidUnderlying LeftComplement(const ArtinBraidUnderlying &b) const
    {
      return Inverse().Product(b);
    };

    ArtinBraidUnderlying RightComplement(const ArtinBraidUnderlying &b) const
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
    // Underlying is the data structure that actually represents the factor (e.g., a permutation table for a braid canonical factor).
    U Underlying;

  public:
    // Factor(under) initializes a new factor, with underlying element under.
    Factor(const U &under)
        : Underlying(under) {}

    // Copy constructor.
    Factor(const Factor &a) const
        : Underlying(a.Underlying) {}

    ~Factor() {}

    // a.OfString sets a to the factor specified by str.
    void OfString(std::string &str)
    {
      Underlying.OfString(str);
    }

    // a.Print(os) prints a to os.
    void Print(std::ostream &os)
    {
      Underlying.Print(os);
    }

    // a.Neutral sets a to Neutral.
    void Neutral()
    {
      Underlying.Neutral();
    }

    // a.Delta() sets a to Delta.
    void Delta()
    {
      Underlying.Delta();
    }

    // a.Compare(b) returns true if a and b are equal, false otherwise.
    bool Compare(const Factor &b) const
    {
      return Underlying.Compare(b.Underlying);
    }

    // a == b returns true if a and b are equal, false otherwise.
    // Syntactic sugar for a.Compare(b).
    bool operator==(const Factor &b) const
    {
      return Compare(b);
    }

    // a != b returns true if a and b are not equal, false otherwise.
    bool operator!=(const Factor &b) const
    {
      return !Compare(b);
    }

    // a.IsDelta() returns whether a = e.
    bool IsNeutral() const
    {
      Factor e = Factor(this);
      e.Neutral();
      return Compare(e);
    }

    // a.IsDelta() returns whether a = Delta.
    bool IsDelta() const
    {
      Factor delta = Factor(this);
      delta.Delta();
      return Compare(delta);
    }

    // a.LeftComplement(b) returns (assuming that a right-divides b) the left complement of a under b, ba^{-1}.
    Factor LeftComplement(const Factor &b) const
    {
      return Factor(Underlying.LeftComplement(b.Underlying));
    }

    // a.LeftComplement() return a's left complement.
    Factor LeftComplement() const
    {
      Factor delta = Factor(this);
      delta.Delta();
      return LeftComplement(delta);
    }

    // a.RightComplement(b) returns (assuming that a left-divides b) the right complement of a under b, a^{-1}b.
    Factor RightComplement(const Factor &b) const
    {
      return Factor(Underlying.RightComplement(b.Underlying));
    }

    // a.RightComplement() return a's right complement.
    Factor RightComplement() const
    {
      Factor delta = Factor(this);
      delta.Delta();
      return RightComplement(delta);
    }

    // ~a return a's right complement.
    // Syntactic sugar for a.RightComplement().
    Factor operator~() const
    {
      return RightComplement();
    }

    // a / b returns (assuming that b left-divides a) the right complement of a under b, b^{-1}a.
    // Syntactic sugar for b.RightComplement(a).
    Factor operator/(const Factor &b) const
    {
      return b.RightComplement(*this);
    }

    // a.DeltaConjugate(k) returns a, conjugated by Delta ^ k.
    // Makes 2 |k| complement calculations.
    Factor DeltaConjugate(const sint16 k) const
    {
      Factor conjugate = Factor(this);
      sint16 i;
      if (k >= 0)
      {
        for (i = 0; i < k; i++)
        {
          conjugate = conjugate.RightComplement().RightComplement();
        }
        return conjugate;
      }
      else
      {
        for (i = 0; i > k; i--)
        {
          conjugate = conjugate.LeftComplement().LeftComplement();
        }
        return conjugate;
      }
    }

    // a.DeltaConjugate() returns a conjugated by Delta.
    // This is done by taking complements twice.
    Factor DeltaConjugate() const
    {
      DeltaConjugate(0);
    }

    // a.LeftMeet(b) returns the left meet of a and b.
    Factor LeftMeet(const Factor &b) const
    {
      return Factor(Underlying.LeftMeet(b.Underlying));
    }

    // a ^ b returns the left meet of a and b.
    // Syntactic sugar for a.LeftMeet(b).
    Factor operator^(const Factor &b) const
    {
      return LeftMeet(b);
    }

    // a.RightMeet(b) returns the right meet of a and b.
    Factor RightMeet(const Factor &b) const
    {
      return Factor(Underlying.RightMeet(b.Underlying));
    }

    // a.LeftJoin(b) returns the left join of a and b.
    Factor LeftJoin(const Factor &b) const {
        return RightComplement().RightMeet(b.RightComplement()).LeftComplement()}

    // a.RightJoin(b) returns the right join of a and b.
    Factor RightJoin(const Factor &b) const
    {
      return LeftComplement().LeftMeet(b.LeftComplement()).RightComplement()
    }

    // a.IsLeftWeighted(b) returns true if a | b is left weighted, or false otherwise.
    bool IsLeftWeighted(const Factor &b) const
    {
      RightComplement().LeftMeet(b).IsTrivial();
    }

    // a.IsRightWeighted(b) returns true if a | b is right weighted, or false otherwise.
    bool IsRightWeighted(const Factor &b) const
    {
      LeftMeet(b.LeftComplement()).IsTrivial();
    }

    // a.Product(b) returns the product of two factors, under the assumption that it lies below Delta.
    Factor Product(const Factor &b) const
    {
      return Factor(Underlying.Product(b.Underlying));
    }

    // a * b is the product of a and b, under the assumption that it lies below Delta.
    // Syntactic sugar for a.Product(b).
    Factor operator*(const Factor &b) const
    {
      return Product(b);
    }

    // a.Randomize() sets a to a random factor.
    void Randomize()
    {
      Underlying.Randomize();
    }

    // a.Atoms() returns the list of the atoms.
    // Specific information about the group (e.g. the number of strands for braids) is recovered by the corresponding U method.
    std::list<Factor> Atoms() const
    {
      std::list<U> atoms = Underlying.Atoms();
      typename std::list<U>::iterator atoms_it;
      std::list<Factor> factor_atoms;
      for (auto const &atoms_it : atoms)
      {
        Factor f = Factor(this);
        f.Set(atoms_it);
        factor_atoms.push_back(f);
      }
      return factor_atoms;
    }
  };

  // MakeLeftWeighted(u, v) computes the left-weighted decomposition u' | v' = u | v, and sets u = u' and v = v'.
  template <class U>
  void MakeLeftWeighted(Factor<U> &u, Factor<U> &v)
  {
    Factor<U> t = ~u ^ v;
    v = v / t;
    u = u * t;
  }

  // MakeRightWeighted(u, v) computes the right-weighted decomposition u' | v' = u | v, and sets u = u' and v = v'.
  template <class U>
  void MakeRightWeighted(Factor<U> &u, Factor<U> &v)
  {
    Factor<U> t = u.RightMeet(v.LeftComplement());
    u = t.LeftComplement(u);
    v = t * v;
  }

  typedef Factor<ArtinBraidUnderlying> ArtinBraidFactor;

  template <class F>
  class Braid
  {

  public:
    // `LeftDelta` is the number of Deltas on the left end of the word.
    sint32 LeftDelta;

    // `RightDelta` is the number of Deltas on the right end of the word.
    sint32 RightDelta;

    // `FactorList` is a list a canonical factors.
    std::list<F> FactorList;

  public:
    // Iterator types for canonical factors.
    typedef typename std::list<F>::iterator FactorItr;
    typedef typename std::list<F>::const_iterator ConstFactorItr;
    typedef typename std::list<F>::reverse_iterator RevFactorItr;
    typedef typename std::list<F>::const_reverse_iterator ConstRevFactorItr;

  public:
    // Constructor that creates an empy word.
    Braid()
        : LeftDelta(0),
          RightDelta(0),
          FactorList() {}

    // Copy constructor.
    Braid(const Braid &w)
        : LeftDelta(w.LeftDelta),
          RightDelta(w.RightDelta),
          FactorList(w.FactorList) {}

    // Constructor that transforms a factor into a word.
    // If the factor is Delta, LeftDelta is incremented.
    Braid(const F &f)
    {
      if f
        .IsDelta()
        {
          LeftDelta = 1;
        }
      else
      {
        FactorList.push_back(f);
      }
    }

    // Destructor.
    ~Braid() {}

    // `w.Neutral()` sets w to the empty word.
    void Neutral()
    {
      LeftDelta = 0;
      RightDelta = 0;
      FactorList.clear();
    }

    // `u.Compare(v)` returns whether u and v have the same internal representation.
    // Should hence only be used on words in normal form.
    bool Compare(Braid &v) const
    {
      return (LeftDelta == v.LeftDelta &&
              RightDelta == v.RightDelta &&
              FactorList == v.FactorList);
    }

    // `u == v` returns whether u and v have the same internal representation.
    // Should hence only be used on words in normal form.
    // Syntactic sugar for `u.Compare(v)`.
    bool operator==(const Braid &v) const
    {
      return Compare(v);
    }

    // `u != v` returns whether u and v do not have the same internal representation.
    // Should hence only be used on words in normal form.
    bool operator!=(const Braid &v) const
    {
      return !Compare(v);
    }

    // `u.IsNeutral` returns whether u represents the neutral element.
    bool IsNeutral() const
    {
      Braid e;
      e.Neutral();
      return Compare(e)
    }

    // `u.Inverse()` returns the inverse of u, not necessarily in a canonical form.
    Braid Inverse() const
    {
      Braid b;
      b.LeftDelta = - RightDelta;
      b.RightDelta = 0;
      for (ConstRevFactorItr it = FactorList.rbegin();
           it != FactorList.rend();
           ++it)
      {
        // Rewrite a_1...a_k Delta^r (*it)^(-1) as
        // a_1...a_k (Delta^r f Delta^(-r)) Delta^(r-1).
        b.FactorList.push_back((*it).RightComplement().DeltaConjugate(-b.RightDelta));
        --b.RightDelta;
      }
      b.RightDelta -= LeftDelta;
      return b;
    }

    // `!u` returns the inverse of u, not necessarily in a canonical form.
    // Syntactic sugar for `u.Inverse()`.
    Braid operator!() const
    {
      return Inverse();
    }

    // `u.LeftProduct(f)` assigns fu, not necessarily in a canonical form, to u.
    void LeftProduct(const F &f) {
      FactorList.push_front(f.DeltaConjugate(LeftDelta));
    }

    // `u.RightProduct(f)` assigns uf, not necessarily in a canonical form, to u.
    void RightProduct(const F &f) {
      FactorList.push_back(f.DeltaConjugate(- RightDelta));
    }

    // `u.Product(v)` returns uv, not necessarily in a canonical form.
    Braid Product(const Braid &v) const {
      Braid w(this);
      w.RightDelta += v.LeftDelta;
      for (ConstFactorItr it = v.FactorList.begin();
           it != v.FactorList.end();
           ++it)
      {
        w.RightMultiply(*it);
      }
      w.RightDelta += v.RightDelta;
      return w;
    }

    // `u * v` returns uv, not necessarily in a canonical form.
    // Syntactic sugar for `u.Product(v)`.
    Braid operator*(const Braid &v) const {
      return Product(v);
    }

  };

}
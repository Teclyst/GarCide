#ifndef CGARSIDE
#define CGARSIDE 

#include <list>
#include <string>
#include <iostream>
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

  static sint16 MadeLeftWeighted = 0;

  // Exceptions.
  struct InvalidStringError
  {
  };
  struct NonMatchingIndexes
  {
  };

  template <class ForItr, class BinFunc>
  inline ForItr apply_binfun(ForItr first, ForItr last, BinFunc f)
  {
    ForItr i, j;
    if ((i = j = first) == last)
      return first;
    while (++j != last && f(*(i++), *j))
      ;
    return j;
  }

  template <class BiItr, class BinFunc>
  inline BiItr reverse_apply_binfun(BiItr first, BiItr last, BinFunc f)
  {
    BiItr i, j;
    if (first == (i = j = last))
      return first;
    --i;
    while ((j = i) != first && f(*--i, *j))
      ;
    return --j;
  }

  template <class ForItr, class BinFun>
  inline void bubble_sort(ForItr first, ForItr last, BinFun f)
  {
    ForItr i;
    if (first == (i = last))
      return;
    while (i != first)
      apply_binfun(--i, last, f);
  }

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
    Factor(const Factor &a)
        : Underlying(a.Underlying) {}

    ~Factor() {}

    // a.OfString sets a to the factor specified by str.
    void OfString(std::string &str)
    {
      Underlying.OfString(str);
    }

    // a.Debug(os) prints a's internal representation to os.
    void Debug(std::ostream &os)
    {
      Underlying.Debug(os);
    }

    // a.Print(os) prints a to os.
    void Print(std::ostream &os)
    {
      Underlying.Print(os);
    }

    // a.Identity sets a to Identity.
    void Identity()
    {
      Underlying.Identity();
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

    Factor &Assign(const Factor &b)
    {
      Underlying.Assign(b.Underlying);
      return *this;
    }

    Factor &operator=(const Factor &b)
    {
      return Assign(b);
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

    // a.IsDelta() returns whether a == e.
    bool IsIdentity() const
    {
      Factor e = Factor(*this);
      e.Identity();
      return Compare(e);
    }

    // a.IsDelta() returns whether a = Delta.
    bool IsDelta() const
    {
      Factor delta = Factor(*this);
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
      Factor delta = Factor(*this);
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
      Factor delta = Factor(*this);
      delta.Delta();
      return RightComplement(delta);
    }

    // ~a return a's right complement.
    // Syntactic sugar for a.RightComplement().
    Factor operator~() const
    {
      return RightComplement();
    }

    // Syntactic sugar for b.RightComplement(a).
    Factor operator/(const Factor &b) const
    {
      return b.RightComplement(*this);
    }

    // a.DeltaConjugate(k) returns a, conjugated by Delta ^ k.
    // Makes 2 |k| complement calculations.
    Factor DeltaConjugate(const sint16 k) const
    {
      Factor conjugate = Factor(*this);
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
    Factor LeftJoin(const Factor &b) const
    {
      return RightComplement().RightMeet(b.RightComplement()).LeftComplement();
    }

    // a.RightJoin(b) returns the right join of a and b.
    Factor RightJoin(const Factor &b) const
    {
      return LeftComplement().LeftMeet(b.LeftComplement()).RightComplement();
    }

    // a.IsLeftWeighted(b) returns true if a | b is left weighted, or false otherwise.
    bool IsLeftWeighted(const Factor &b) const
    {
      return RightComplement().LeftMeet(b).IsTrivial();
    }

    // a.IsRightWeighted(b) returns true if a | b is right weighted, or false otherwise.
    bool IsRightWeighted(const Factor &b) const
    {
      return LeftMeet(b.LeftComplement()).IsTrivial();
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
        Factor f = Factor(*this);
        f.Set(atoms_it);
        factor_atoms.push_back(f);
      }
      return factor_atoms;
    }
  };

  // MakeLeftWeighted(u, v) computes the left-weighted decomposition u' | v' = u | v, and sets u = u' and v = v'.
  // It then returns true if something was done (so that it may be used with `apply_binfun`).
  template <class F>
  bool MakeLeftWeighted(F &u, F &v)
  {
    MadeLeftWeighted++;
    F t = (~u) ^ v;
    if (t.IsIdentity())
    {
      return false;
    }
    else
    {
      v = v / t;
      u = u * t;
      return true;
    }
  }

  // We maintain braids in LCF at all time.
  template <class F>
  class Braid
  {

  public:
    // `Delta` is the number of Deltas on the left end of the word.
    sint32 Delta;

    // `FactorList` is a list a canonical factors.
    std::list<F> FactorList;

  public:
    // Iterator types for canonical factors.
    typedef typename std::list<F>::iterator FactorItr;
    typedef typename std::list<F>::const_iterator ConstFactorItr;
    typedef typename std::list<F>::reverse_iterator RevFactorItr;
    typedef typename std::list<F>::const_reverse_iterator ConstRevFactorItr;

  public:
    // Constructor that creates an empty word.
    Braid()
        : Delta(0),
          FactorList() {}

    // Constructor, with expected canonical length.
    Braid(const sint16 expected_length)
        : Delta(0),
          FactorList(expected_length) {}

    // Copy constructor.
    Braid(const Braid &w)
        : Delta(w.Delta),
          FactorList(w.FactorList) {}

    // Constructor that transforms a factor into a word.
    // If the factor is Delta, `Delta` is incremented.
    Braid(const F &f)
    {
      if (f
              .IsDelta())
      {
        Delta = 1;
      }
      else
      {
        FactorList.push_back(f);
      }
    }

    // Destructor.
    ~Braid() {}

    void Debug(std::ostream &os) {
      os << Delta;
      for(FactorItr it = FactorList.begin(); it != FactorList.end(); it++){
        (*it).Debug(os);
        os << " ";
      }
    }

    // `w.Identity()` sets w to the empty word.
    void Identity()
    {
      Delta = 0;
      FactorList.clear();
    }

    // `u.Compare(v)` returns whether u and v have the same internal representation.
    bool Compare(Braid &v) const
    {
      return (Delta == v.Delta &&
              FactorList == v.FactorList);
    }

    // `u == v` returns whether u and v have the same internal representation.
    // Syntactic sugar for `u.Compare(v)`.
    bool operator==(const Braid &v) const
    {
      return Compare(v);
    }

    // `u != v` returns whether u and v do not have the same internal representation.
    bool operator!=(const Braid &v) const
    {
      return !Compare(v);
    }

    // `u.IsIdentity` returns whether u represents the Identity element.
    bool IsIdentity() const
    {
      return Compare(Delta == 0 && FactorList.empty());
    }

    // `u.Inverse()` returns the inverse of u, not necessarily in a canonical form.
    Braid Inverse() const
    {
      Braid b;
      b.Delta = Delta;
      for (ConstFactorItr it = FactorList.begin();
           it != FactorList.end();
           it++)
      {
        // Rewrite a_1...a_k (f)^(- 1) Delta^r as
        // a_1...a_k Delta ^ (r - 1) (Delta^(- r) d_L(f) Delta^r).
        b.FactorList.push_front((*it).RightComplement().DeltaConjugate(b.Delta));
        --b.Delta;
      }
      return b;
    }

    // `!u` returns the inverse of u, not necessarily in a canonical form.
    // Syntactic sugar for `u.Inverse()`.
    Braid operator!() const
    {
      return Inverse();
    }

    // Clean gets rid of (factor) Deltas at the start, and identity elements at the end pf `FactorList`.
    void Clean()
    {
      FactorItr it = FactorList.begin();
      while (it != FactorList.end() && (*it).IsDelta())
      {
        ++it;
        ++Delta;
      }
      FactorList.erase(FactorList.begin(), it);
      RevFactorItr revit = FactorList.rbegin();
      while (revit != FactorList.rend() && (*revit).IsIdentity())
      {
        ++revit;
      }
      FactorList.erase(revit.base(), FactorList.end());
    }

    // `u.LeftProduct(f)` assigns fu to u.
    void LeftProduct(F &f)
    {
      FactorList.push_front(f.DeltaConjugate(Delta));
      apply_binfun(FactorList.begin(), FactorList.end(), MakeLeftWeighted<F>);
      Clean();
    }

    // `u.RightProduct(f)` assigns uf, not necessarily in a canonical form, to u.
    void RightProduct(F &f)
    {
      FactorList.push_back(f);
      reverse_apply_binfun(FactorList.begin(), FactorList.end(), MakeLeftWeighted<F>);
      Clean();
    }

    // `u.Product(v)` returns uv, not necessarily in a canonical form.
    Braid Product(const Braid &v) const
    {
      Braid w(this);
      w.Delta += v.Delta;
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
    Braid operator*(const Braid &v) const
    {
      return Product(v);
    }

    void Normalize()
    {
      bubble_sort(FactorList.begin(), FactorList.end(), MakeLeftWeighted<F>);
      Clean();
    }

    // Randomizes the braid, setting it at a given length, using parameters extracted from Factor f.
    // The result isn't in LCF.
    void Randomize(const F &f, sint16 canonical_length)
    {
      FactorList.clear(); // Potential memory leak? To be checked.
      Delta = 0;
      for (sint16 i = 0; i < canonical_length; i++)
      {
        F g = F(f);
        g.Randomize();
        FactorList.push_back(g);
      }
    }

    void Print(std::ostream &os)
    {
      os << "D^" << Delta << " . ";
      for (FactorItr it = FactorList.begin(); it != FactorList.end(); it++)
      {
        (*it).Print(os);
        os << " . ";
      }
      os << std::endl;
    }
  };

}

#endif
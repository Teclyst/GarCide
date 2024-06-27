#include "cgarside.h"

namespace SSS
{
  using namespace CGarside;

  template <class F>
  Braid<F> SendToSSS(const Braid<F> b)
  {

    typename F::ParameterType n = b.GetParameter();

    sint16 k = F(n).LatticeHeight();

    Braid<F> b2 = Braid(b), b3 = Braid(b);

    sint16 p = b.Delta;
    sint16 j = 0;

    while (j <= k)
    {
      b2.Cycling();

      if (b2.LeftDelta == p)
      {
        j++;
      }
      else
      {
        b3 = b2;
        p++;
        j = 0;
      }
    }

    j = 0;
    sint16 l = b3.Sup();
    while (j <= k)
    {
      b3.Decycling();

      if (b3.Sup() == l)
      {
        j++;
      }
      else
      {
        b2 = b3;
        l--;
        j = 0;
      }
    }

    return b2;
  }

  template <class F>
  Braid<F> SendToSSS(const Braid<F> b, Braid<F> &c)
  {

    typename F::ParameterType n = b.GetParameter();

    sint16 k = F(n).LatticeHeight();

    Braid<F> b2 = Braid(b), b3 = Braid(b), c2 = Braid(n);

    c.Identity();

    sint16 p = b.Delta;
    sint16 j = 0;

    while (j <= k)
    {
      if (b2.CanonicalLength() == 0)
      {
        return b2;
      }

      c2.RightProduct(b2.FactorList.front().DeltaConjugate(b2.delta));
      b2.Cycling();

      if (b2.LeftDelta == p)
      {
        j++;
      }
      else
      {
        b3 = b2;
        p++;
        j = 0;
        c.RightProduct(c2);
        c2.Identity();
      }
    }

    j = 0;
    sint16 l = b3.Sup();
    c2.Identity();

    while (j <= k)
    {
      c2.LeftProduct(b3.FactorList.back());
      b3.Decycling();

      if (b3.Sup() == l)
      {
        j++;
      }
      else
      {
        b2 = b3;
        l--;
        j = 0;
        c.RightDivide(c2);
        c2.Identity();
      }
    }

    return b2;
  }

  template <class F>
  F MinSS(const Braid<F> b, const F f)
  {
    F r2 = F(f), r = F(f);
    r.Identity();

    Braid<F> w = Braid(b);
    w.Delta = 0;

    while (!r2.IsIdentity())
    {
      r.RightProduct(r2);
      r2 = (w * r).Remainder(R.DeltaConjugate(b.Delta))
    }

    return r;
  }

  template <class F>
  F MinSSS(const Braid<F> b, const F f)
  {
    F r = MinSS(b, f);
    Braid<F> b2 = Braid(b);
    b2.Conjugate(r);
    b2.MakeRCFFromLCF();
    while (b2.CanonicalLength() > b.CanonicalLength())
    {
      r.RightProduct(b2.FactorList.front());
      b2.ConjugateRCF(r);
    }
    return r;
  }

  template <class F>
  std::vector<F> MinSSS(Braid<F> b)
  {
    F f = F(b.GetParameter());
    std::list<F> atoms = f.Atoms();

    std::vector<F> min;

    static bool table[CGarside::MaxBraidIndex];
    bool should_be_added;

    for (sint16 i = 0; i < atoms.size(); i++)
    {
      f = MinSSS(b, atoms[i]);
      should_be_added = true;

      // We check, before adding f, that a divisor of it wasn't added already with some other atom dividing it.
      for (sint16 j = 0; j < i && should_be_added; j++)
      {
        should_be_added = !(table[j] && atoms[j] ^ f == atoms[j]);
      }
      // If that is not the case, we also check if the atom we see is the last that might generate f.
      // This is to avoid duplicates; furthermore, if some strict left divisor of f can be generated by an atom,
      // doing things in this order ensures that it would have been detected by the first loop by the time we try to add f.
      for (sint16 j = i + 1; j < atoms.size() && should_be_added; j++)
      {
        should_be_added = !(atoms[j] ^ f == atoms[j]);
      }
      if (should_be_added)
      {
        min.push_back(f);
        table[i] = true;
      }
    }

    return min;
  }

  template <class F>
  std::unordered_set<Braid<F>> SSS(Braid<F> b)
  {
    Braid<F> b2 = SendToSSS(b);
    F f = F(b.GetParameter());

    std::list<Braid<F>> queue;
    std::unordered_set<Braid<F>> sss;

    queue.push_back(Braid(b2));
    sss.insert(b2);

    while (!queue.empty())
    {
      Braid b3 = queue.front();
      queue.pop_front();

      std::vector<F> min = MinSSS(b3);

      for (typename std::vector<F>::iterator it = min.begin();
           it != min.end();
           it++)
      {
        b2 = b3;
        b2.Conjugate(*it);
        if (sss.find(b2) == sss.end())
        {
          sss.insert(b2);
          queue.push_back(Braid(b2));
        }
      }
    }

    return sss;
  }

  template <class F>
  inline bool AreConjugate(Braid<F> u, Braid<F> v)
  {
    std::unordered_set<Braid<F>> u_sss = SSS(u);
    return (u_sss.find(SendToSSS(v)) != u_sss.end());
  }

}

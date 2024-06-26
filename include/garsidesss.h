#include "cgarside.h"
#include <unordered_set>

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
  std::unordered_set<F> MinSSS(Braid<F> b)
  {
    sint16 i, j, k;

    // foo is by itself useless, and only built because Atoms is a method.
    F foo = F(b.GetParameter());
    std::list<F> atoms = foo.Atoms();

    std::unordered_set<F> min;

    for (typename Braid<F>::FactorItr it = atoms.begin(); it != atoms.end; it++)
    {
      min.insert(MinSSS(b, *it));
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

      std::unordered_set<F> min = MinSSS(b3);

      for (typename std::unordered_set<F>::iterator it = min.begin();
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

  template<class F>
  inline bool AreConjugate(Braid<F> u, Braid<F> v) {
    std::unordered_set<Braid<F>> u_sss = SSS(u);
    return (u_sss.find(SendToSSS(v)) != u_sss.end());
  }

}

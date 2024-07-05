/*
    Copyright (C) 2004 Juan Gonzalez-Meneses.

    This file is part of Braiding.

    Braiding is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    Braiding is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Braiding; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

/*
    braiding.cpp,  v 1.0.   04/10/2004
    Juan Gonzalez-Meneses <meneses(at)us.es>
*/

#include "cbraid.h"
#include <iostream>
#include <iomanip>
#include <fstream>

namespace Braiding
{

  using CBraid::ArtinBraid;
  using CBraid::ArtinFactor;
  using CBraid::sint16;
  using std::cout;
  using std::endl;
  using std::list;

  // typedef ArtinPresentation P;

  ///////////////////////////////////////////////////////
  //
  //  CL(B)  computes the Canonical length of a braid B,
  //         given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  sint16 CL(CBraid::Braid<P> B)
  {
    sint16 n = 0;
    typename CBraid::Braid<P>::ConstFactorItr it;

    for (it = B.FactorList.begin(); it != B.FactorList.end(); it++)
      n++;

    return n;
  };

  ///////////////////////////////////////////////////////
  //
  //  Sup(B)  computes the supremun of a braid B,
  //          given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  sint16 Sup(CBraid::Braid<P> B)
  {
    sint16 s;

    s = CL(B) + B.LeftDelta;

    return s;
  }

  ///////////////////////////////////////////////////////
  //
  //  Cycling(B)  computes the cycling of a braid B,
  //              given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> Cycling(CBraid::Braid<P> B)
  {

    sint16 n;

    if (CL(B) == 0)
      return B;

    //   B.MakeLCF();

    n = B.Index();

    CBraid::Factor<P> F = CBraid::Factor<P>(n);

    F = *B.FactorList.begin();
    B.FactorList.push_back(F.Flip(-B.LeftDelta));
    B.FactorList.pop_front();
    B.MakeLCF();

    return B;
  }

  ///////////////////////////////////////////////////////
  //
  //  Decycling(B)  computes the decycling of a braid B,
  //                given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> Decycling(CBraid::Braid<P> B)
  {
    sint16 n;

    if (CL(B) == 0)
      return B;

    //   B.MakeLCF();

    n = B.Index();

    CBraid::Factor<P> F = CBraid::Factor<P>(n);

    F = B.FactorList.back();
    B.FactorList.push_front(F.Flip(B.LeftDelta));
    B.FactorList.pop_back();
    B.MakeLCF();

    return B;
  }

  /////////////////////////////////////////////////////////////
  //
  //  WordToBraid(w,n)  Transforms a word w (list of letters)
  //                    into a braid on n strands in LCF.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> WordToBraid(list<sint16> w, sint16 n)
  {
    CBraid::Braid<P> B = CBraid::Braid<P>(n);
    CBraid::Braid<P> B2 = CBraid::Braid<P>(n);
    CBraid::Factor<P> F = CBraid::Factor<P>(n);
    sint16 k, sigma, i;
    list<sint16>::iterator it;

    for (it = w.begin(); it != w.end(); it++)
    {
      if (*it > 0)
      {
        if (*it == n)
          F.Delta();
        else
        {
          F.Identity();
          sigma = *it;
          k = F[sigma];
          F[sigma] = F[sigma + 1];
          F[sigma + 1] = k;
        }
        B.RightMultiply(F);
      }
      else
      {
        if (*it == -n)
        {
          B2.Identity();
          B2.LeftDelta = -1;
        }
        else
        {
          F.Identity();
          i = -*it;
          k = F[i];
          F[i] = F[i + 1];
          F[i + 1] = k;
          F = (~F).Flip();
          B2 = (!CBraid::Braid<P>(CBraid::Factor<P>(n, 1))) * CBraid::Braid<P>(F);
        }
        B.RightMultiply(B2);
      }
    }
    B.MakeLCF();
    return B;
  }

  /////////////////////////////////////////////////////////////
  //
  //  PrintBraidWord(B)  Shows on the screen the braid B
  //                     (given in left or right normal form)
  //                     written as a word in Artin generators.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  void PrintBraidWord(CBraid::Braid<P> B)
  {
    if (B.LeftDelta == 1)
    {
      cout << "D";
      if (CL(B))
        cout << " . ";
    }
    else if (B.LeftDelta != 0)
    {
      cout << "D^(" << B.LeftDelta << ")";
      if (CL(B))
        cout << " . ";
    }
    sint16 i, j, k, n = B.Index();
    CBraid::Factor<P> F = CBraid::Factor<P>(n);

    typename list<CBraid::Factor<P>>::iterator it;

    for (it = B.FactorList.begin(); it != B.FactorList.end(); it++)
    {
      if (it != B.FactorList.begin())
        cout << ". ";

      (*it).Print(cout);
    }

    if (B.RightDelta == 1)
    {
      if (CL(B))
        cout << ". ";
      cout << "D";
    }
    else if (B.RightDelta != 0)
    {
      if (CL(B))
        cout << ". ";
      cout << "D^(" << B.RightDelta << ")";
    }
  }

  /////////////////////////////////////////////////////////////
  //
  //  PrintBraidWord(B,file)  Prints on 'file' the braid B (given in LCF)
  //                          written as a word in Artin generaotrs.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  void PrintBraidWord(CBraid::Braid<P> B, char *file)
  {
    std::ofstream f(file, std::ios::app);

    if (B.LeftDelta == 1)
    {
      f << "D";
      if (CL(B))
        f << " . ";
    }
    else if (B.LeftDelta != 0)
    {
      f << "D^(" << B.LeftDelta << ")";
      if (CL(B))
        f << " . ";
    }
    sint16 i, j, k, n = B.Index();
    CBraid::Factor<P> F = CBraid::Factor<P>(n);

    typename list<CBraid::Factor<P>>::iterator it;

    for (it = B.FactorList.begin(); it != B.FactorList.end(); it++)
    {
      if (it != B.FactorList.begin())
        f << ". ";

      (*it).Print(f);
    }

    if (B.RightDelta == 1)
    {
      if (CL(B))
        f << ". ";
      f << "D";
    }
    else if (B.RightDelta != 0)
    {
      if (CL(B))
        f << ". ";
      f << "D^(" << B.RightDelta << ")";
    }
    f.close();
  }

  /////////////////////////////////////////////////////////////
  //
  //  PrintWord(word,n,power)  Shows on the screen the braid "word"
  //                           on n strands raised to some "power".
  //
  /////////////////////////////////////////////////////////////

  void PrintWord(list<sint16> &word, sint16 n, sint16 power)
  {
    list<sint16>::iterator itw;

    if (power != 1)
      cout << "( ";

    for (itw = word.begin(); itw != word.end(); itw++)
    {
      if (*itw == n)
        cout << "D ";
      else if (*itw == -n)
        cout << "-D ";
      else
        cout << *itw << " ";
    }

    if (power != 1)
      cout << ")^" << power;
  }

  /////////////////////////////////////////////////////////////
  //
  //  PrintWord(word,n,power,file)  Prints on "file" the braid "word"
  //                                on n strands raised to some "power".
  //
  /////////////////////////////////////////////////////////////

  void PrintWord(list<sint16> &word, sint16 n, sint16 power, char *file)
  {
    list<sint16>::iterator itw;
    std::ofstream f(file, std::ios::app);

    if (power != 1)
      f << "( ";

    for (itw = word.begin(); itw != word.end(); itw++)
    {
      if (*itw == n)
        f << "D ";
      else if (*itw == -n)
        f << "-D ";
      else
        f << *itw << " ";
    }

    if (power != 1)
      f << ")^" << power;

    f.close();
  }

  /////////////////////////////////////////////////////////////
  //
  //  Crossing(word,n,power)  Computes the crossing numbers of
  //                           the braid on n strands given by
  //				     "word" raised to "power".
  //
  /////////////////////////////////////////////////////////////

  void Crossing(list<sint16> word, sint16 n, sint16 power, sint16 **cross)
  {
    sint16 i, j, k, l, m;
    list<sint16>::iterator itw;
    sint16 *perm = new sint16[n];
    for (i = 1; i <= n; i++)
      perm[i] = i;

    for (i = 1; i < n; i++)
    {
      for (j = i + 1; j <= n; j++)
        cross[i][j] = 0;
    }

    for (m = 1; m <= power; m++)
    {
      for (itw = word.begin(); itw != word.end(); itw++)
      {
        if (*itw == n || *itw == -n)
        {
          if (*itw == n)
            l = 1;
          else
            l = -1;
          for (i = 1; i < n; i++)
          {
            for (j = i + 1; j <= n; j++)
              cross[i][j] += l;
          }
          for (i = 1; 2 * i <= n; i++)
          {
            k = perm[i];
            perm[i] = perm[n + 1 - i];
            perm[n + 1 - i] = k;
          }
        }
        else
        {
          if (*itw > 0)
            l = *itw;
          else
            l = -(*itw);
          if (perm[l] < perm[l + 1])
          {
            i = perm[l];
            j = perm[l + 1];
          }
          else
          {
            i = perm[l + 1];
            j = perm[l];
          }
          if (*itw > 0)
            cross[i][j]++;
          else
            cross[i][j]--;
          k = perm[l];
          perm[l] = perm[l + 1];
          perm[l + 1] = k;
        }
      }
    }
    delete[] perm;
  }

  /////////////////////////////////////////////////////////////
  //
  //  SendToSSS(B)  Computes a braid conjugate to B that
  //                belongs to its Super Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> SendToSSS(CBraid::Braid<P> B)
  {
    sint16 n, k, j, p, l;

    n = B.Index();
    k = n * (n - 1) / 2;

    CBraid::Braid<P> B2 = CBraid::Braid<P>(n), B3 = CBraid::Braid<P>(n);

    B.MakeLCF();

    j = 0;
    B2 = B;
    B3 = B;
    p = B.LeftDelta;

    while (j <= k)
    {
      B2 = Cycling(B2);

      if (B2.LeftDelta == p)
        j++;
      else
      {
        B3 = B2;
        p++;
        j = 0;
      }
    }

    j = 0;
    B2 = B3;
    l = Sup(B2);
    while (j <= k)
    {
      B2 = Decycling(B2);

      if (Sup(B2) == l)
        j++;
      else
      {
        B3 = B2;
        l--;
        j = 0;
      }
    }

    return B3;
  }

  /////////////////////////////////////////////////////////////
  //
  //  SendToSSS(B,C)  Computes a braid conjugate to B that
  //                  belongs to its Super Summit Set, and a braid
  //                  C that conjugates B to the result.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> SendToSSS(CBraid::Braid<P> B, CBraid::Braid<P> &C)
  {
    sint16 n, k, j, p, l;

    n = B.Index();
    k = n * (n - 1) / 2;

    CBraid::Braid<P> B2 = CBraid::Braid<P>(n), B3 = CBraid::Braid<P>(n), C2 = CBraid::Braid<P>(n);

    B.MakeLCF();
    C = CBraid::Braid<P>(n);

    j = 0;
    B2 = B;
    B3 = B;
    p = B.LeftDelta;

    while (j <= k)
    {
      if (CL(B2) == 0)
      {
        C.MakeLCF();
        return B2;
      }
      C2 = C2 * ((*B2.FactorList.begin()).Flip(B2.LeftDelta));
      B2 = Cycling(B2);

      if (B2.LeftDelta == p)
        j++;
      else
      {
        B3 = B2;
        p++;
        j = 0;
        C = C * C2;
        C2 = CBraid::Braid<P>(n);
      }
    }

    j = 0;
    B2 = B3;
    l = Sup(B2);
    C2 = CBraid::Braid<P>(n);
    while (j <= k)
    {
      C2.LeftMultiply(B2.FactorList.back());
      B2 = Decycling(B2);

      if (Sup(B2) == l)
        j++;
      else
      {
        B3 = B2;
        l--;
        j = 0;
        C = C * (!C2);
        C2 = CBraid::Braid<P>(n);
      }
    }

    C.MakeLCF();
    return B3;
  }

  /////////////////////////////////////////////////////////////
  //
  //  LeftWedge(F1,F2)  Given two simple factors F1 and F2, computes
  //                    their left lcm. That is, the smallest simple factor
  //                    F such that F1<F and F2<F.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> LeftWedge(CBraid::Factor<P> F1, CBraid::Factor<P> F2)
  {
    return (~RightMeet(~F1, ~F2)).Flip();
  }

  /////////////////////////////////////////////////////////////
  //
  //  RightWedge(F1,F2)  Given two simple factors F1 and F2, computes
  //                    their right lcm. That is, the smallest simple factor
  //                    F such that F>F1 and F>F2.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> RightWedge(CBraid::Factor<P> F1, CBraid::Factor<P> F2)
  {
    return !LeftWedge(!F1, !F2);
  }

  /////////////////////////////////////////////////////////////
  //
  //  Remainder(B,F)   Given a positive braid B in LCF and a simple
  //                   factor F, computes the simple factor S such
  //                   that BS=LeftWedge(B,F).
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> Remainder(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Factor<P> Fi = F;
    if (B.LeftDelta != 0)
    {
      Fi.Identity();
      return Fi;
    }

    typename list<CBraid::Factor<P>>::iterator it;
    for (it = B.FactorList.begin(); it != B.FactorList.end(); it++)
    {

      Fi = !(*it) * LeftWedge(*it, Fi);
    }
    return Fi;
  }

  /////////////////////////////////////////////////////////////
  //
  //  LeftMeet(B1,B2)  Given two braids B1 and B2, computes
  //                    their left gcd. That is, the greatest braid
  //                    B such that B<B1 and B<B2.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> LeftMeet(CBraid::Braid<P> B1, CBraid::Braid<P> B2)
  {
    sint16 n = B1.Index(), shift = 0;
    CBraid::Braid<P> B = CBraid::Braid<P>(n);
    CBraid::Factor<P> F1 = CBraid::Factor<P>(n, 0), F2 = CBraid::Factor<P>(n, 0), F = CBraid::Factor<P>(n, 1);

    B1.MakeLCF();
    B2.MakeLCF();

    shift -= B1.LeftDelta;
    B2.LeftDelta -= B1.LeftDelta;
    B1.LeftDelta = 0;

    if (B2.LeftDelta < 0)
    {
      shift -= B2.LeftDelta;
      B1.LeftDelta -= B2.LeftDelta;
      B2.LeftDelta = 0;
    }

    while (!F.CompareWithIdentity())
    {
      if (B1.LeftDelta > 0)
        F1 = CBraid::Factor<P>(n, 1);
      else if (CL(B1) == 0)
        F1 = CBraid::Factor<P>(n, 0);
      else
        F1 = *B1.FactorList.begin();

      if (B2.LeftDelta > 0)
        F2 = CBraid::Factor<P>(n, 1);
      else if (CL(B2) == 0)
        F2 = CBraid::Factor<P>(n, 0);
      else
        F2 = *B2.FactorList.begin();

      F = LeftMeet(F1, F2);

      B.RightMultiply(F);
      B1.LeftMultiply(!(CBraid::Braid<P>(F)));
      B1.MakeLCF();
      B2.LeftMultiply(!(CBraid::Braid<P>(F)));
      B2.MakeLCF();
    }

    B.MakeLCF();
    B.LeftDelta -= shift;
    return B;
  }

  /////////////////////////////////////////////////////////////
  //
  //  LeftWedge(B1,B2)  Given two braids B1 and B2, computes
  //                    their left lcm. That is, the smallest braid
  //                    B such that B1<B and B2<B.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> LeftWedge(CBraid::Braid<P> B1, CBraid::Braid<P> B2)
  {
    sint16 n = B1.Index(), shift = 0;
    CBraid::Braid<P> B = CBraid::Braid<P>(n);
    CBraid::Factor<P> F2 = CBraid::Factor<P>(n, 0), F = CBraid::Factor<P>(n, 1);

    B1.MakeLCF();
    B2.MakeLCF();

    shift -= B1.LeftDelta;
    B2.LeftDelta -= B1.LeftDelta;
    B1.LeftDelta = 0;

    if (B2.LeftDelta < 0)
    {
      shift -= B2.LeftDelta;
      B1.LeftDelta -= B2.LeftDelta;
      B2.LeftDelta = 0;
    }

    B = B1;

    while (!B2.CompareWithIdentity())
    {
      if (B2.LeftDelta > 0)
        F2 = CBraid::Factor<P>(n, 1);
      else if (CL(B2) == 0)
        F2 = CBraid::Factor<P>(n, 0);
      else
        F2 = *B2.FactorList.begin();

      F = Remainder(B1, F2);

      B.RightMultiply(F);
      B1.RightMultiply(F);
      B1.LeftMultiply(!(CBraid::Braid<P>(F2)));
      B1.MakeLCF();
      B2.LeftMultiply(!(CBraid::Braid<P>(F2)));
      B2.MakeLCF();
    }

    B.MakeLCF();
    B.LeftDelta -= shift;
    return B;
  }

  /////////////////////////////////////////////////////////////
  //
  //  MinSS(B,F)   Given a braid B in its Summit Set (and in LCF),
  //               computes the minimal simple factor R such that
  //               F<R and B^R is in the Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> MinSS(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Factor<P> R2 = F;
    CBraid::Braid<P> W = B;
    W.LeftDelta = 0;
    CBraid::Factor<P> R = CBraid::Factor<P>(F.Index(), 0);

    while (R2.CompareWithIdentity() == 0)
    {
      R = R * R2;
      R2 = Remainder(W * R, R.Flip(B.LeftDelta));
    }
    return R;
  }

  /////////////////////////////////////////////////////////////
  //
  //  MinSSS(B,F)  Given a braid B in its Super Summit Set (and in LCF),
  //               computes the minimal simple factor R such that
  //               F<R and B^R is in the Super Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> MinSSS(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Factor<P> R = MinSS(B, F);
    sint16 cl = CL(B);
    CBraid::Braid<P> B2 = !(CBraid::Braid<P>(R)) * B * R;
    B2.MakeRCF();
    while (CL(B2) > cl)
    {
      R = R * (*B2.FactorList.begin());
      B2 = (!(CBraid::Braid<P>(R)) * B * R).MakeRCF();
    }
    return R;
  }

  /////////////////////////////////////////////////////////////
  //
  //  MinSSS(B)  Given a braid B in its Super Summit Set (and in LCF),
  //             computes the set of minimal simple factors R that
  //             B^R is in the Super Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Factor<P>> MinSSS(CBraid::Braid<P> B)
  {
    sint16 i, j, k, test;
    sint16 n = B.Index();
    sint16 *table = new sint16[n];

    list<CBraid::Factor<P>> Min;

    for (i = 0; i < n; i++)
      table[i] = 0;
    CBraid::Factor<P> F = CBraid::Factor<P>(n);
    for (i = 1; i < n; i++)
    {
      test = 1;
      F.Identity();
      k = F[i];
      F[i] = F[i + 1];
      F[i + 1] = k;
      F = MinSSS(B, F);
      for (j = 1; j < i; j++)
      {
        if (table[j - 1] && F[j] > F[j + 1])
          test = 0;
      }
      for (j = i + 1; j < n; j++)
      {
        if (F[j] > F[j + 1])
          test = 0;
      }
      if (test)
        Min.push_back(F);
    }
    return Min;
  }

  /////////////////////////////////////////////////////////////
  //
  //  SSS(B)  Given a braid B, computes its Super Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Braid<P>> SSS(CBraid::Braid<P> B)
  {
    CBraid::Braid<P> B2 = SendToSSS(B);
    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index());
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itf;

    list<CBraid::Braid<P>> sss;
    sss.push_back(B2);

    typename list<CBraid::Braid<P>>::iterator it = sss.begin();
    while (it != sss.end())
    {
      Min = MinSSS(*it);
      for (itf = Min.begin(); itf != Min.end(); itf++)
      {
        F = *itf;
        B2 = ((!CBraid::Braid<P>(F)) * (*it) * F).MakeLCF();
        if (find(sss.begin(), sss.end(), B2) == sss.end())
          sss.push_back(B2);
      }
      it++;
    }
    return sss;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Trajectory(B)  Computes the trajectory of a braid B, that is,
  //                 a list containing the iterated cyclings of B,
  //                 until the first repetition.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Braid<P>> Trajectory(CBraid::Braid<P> B)
  {
    list<CBraid::Braid<P>> p;
    typename list<CBraid::Braid<P>>::iterator it;

    while (find(p.begin(), p.end(), B) == p.end())
    {
      p.push_back(B);
      B = Cycling(B);
    }

    return p;
  }

  /////////////////////////////////////////////////////////////
  //
  //  SendToUSS(B)  Computes a braid conjugate to B that
  //                belongs to its Ultra Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> SendToUSS(CBraid::Braid<P> B)
  {
    CBraid::Braid<P> B2 = SendToSSS(B);
    list<CBraid::Braid<P>> T = Trajectory(B2);
    return Cycling(T.back());
  }

  /////////////////////////////////////////////////////////////
  //
  //  SendToUSS(B,C)  Computes a braid conjugate to B that
  //                  belongs to its Ultra Summit Set, and a braid
  //                  C that conjugates B to the result.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> SendToUSS(CBraid::Braid<P> B, CBraid::Braid<P> &C)
  {
    CBraid::Braid<P> B2 = SendToSSS(B, C);
    list<CBraid::Braid<P>> T = Trajectory(B2);

    CBraid::Braid<P> D = Cycling(T.back());

    typename list<CBraid::Braid<P>>::iterator it = T.begin();

    while ((*it) != D)
    {
      C = C * ((*(*it).FactorList.begin()).Flip(B2.LeftDelta));
      it++;
    }

    return D;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Transport(B,F)   Given a braid B (in its USS and in LCF), and a simple factor
  //                   F such that B^F is in its SSS, computes the transport of F.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> Transport(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Braid<P> B2 = ((!CBraid::Braid<P>(F)) * B * F).MakeLCF();
    CBraid::Braid<P> B3 = ((!CBraid::Braid<P>(*B.FactorList.begin())) * F * (*B2.FactorList.begin())).MakeLCF();
    return *B3.FactorList.begin();
  }

  /////////////////////////////////////////////////////////////
  //
  //  Returns(B,F)   Given a braid B (in its USS and in LCF), and a simple factor
  //                 F such that B^F is in its SSS, computes the iterated transports
  //                 of F that send B to an element in the trajectory of B^F, until
  //                 the first repetition.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Factor<P>> Returns(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    list<CBraid::Factor<P>> ret;
    typename list<CBraid::Factor<P>>::iterator it = ret.end();
    sint16 n = B.Index();
    CBraid::Braid<P> B1 = B, C1 = CBraid::Braid<P>(n), C2 = CBraid::Braid<P>(n);
    sint16 i, N = 1;
    CBraid::Factor<P> F1 = F;

    C1 = CBraid::Braid<P>((*B1.FactorList.begin()).Flip(B1.LeftDelta));
    B1 = Cycling(B1);
    while (B1 != B)
    {
      C1.RightMultiply((*B1.FactorList.begin()).Flip(B1.LeftDelta));
      B1 = Cycling(B1);
      N++;
    }

    while (it == ret.end())
    {
      ret.push_back(F1);

      B1 = ((!CBraid::Braid<P>(F1)) * B * F1).MakeLCF();

      C2.Identity();

      for (i = 0; i < N; i++)
      {
        C2.RightMultiply((*B1.FactorList.begin()).Flip(B1.LeftDelta));
        B1 = Cycling(B1);
      }

      CBraid::Braid<P> B2 = ((!C1) * F1 * C2).MakeLCF();

      if (B2.LeftDelta == 1)
        F1 = CBraid::Factor<P>(n, 1);
      else if (B2.CompareWithIdentity())
        F1 = CBraid::Factor<P>(n, 0);
      else
        F1 = *B2.FactorList.begin();

      it = find(ret.begin(), ret.end(), F1);
    }

    while (it != ret.begin())
      ret.pop_front();

    return ret;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Pullback(B,F)   Given a braid B (in its USS and in LCF), and a
  //                  simple factor F, computes the pullback of F.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> Pullback(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Factor<P> F1 = (*B.FactorList.begin());
    F1 = F1.Flip(B.LeftDelta + 1);
    CBraid::Factor<P> F2 = F;
    F2 = F2.Flip();

    CBraid::Braid<P> B2 = CBraid::Braid<P>(F1) * F2;

    CBraid::Factor<P> delta = CBraid::Factor<P>(B.Index(), 1);
    B2 = (B2 * Remainder(B2, delta)).MakeLCF();

    (B2.LeftDelta)--;

    CBraid::Factor<P> b0 = CBraid::Factor<P>(B.Index());

    if (B2.LeftDelta == 1)
      b0 = CBraid::Factor<P>(B.Index(), 1);
    else if (B2.CompareWithIdentity())
      b0 = CBraid::Factor<P>(B.Index(), 0);
    else
      b0 = *B2.FactorList.begin();

    CBraid::Factor<P> bi = F.Flip(B.LeftDelta);

    typename list<CBraid::Factor<P>>::iterator it;
    for (it = B.FactorList.begin(); it != B.FactorList.end(); it++)
    {
      if (it != B.FactorList.begin())
        bi = (!(*it)) * LeftWedge(bi, *it);
    }
    CBraid::Factor<P> b = LeftWedge(b0, bi);
    return MinSSS(B, b);
  }

  /////////////////////////////////////////////////////////////
  //
  //  MainPullback(B,F)   Given a braid B (in its USS and in LCF), and a
  //                      simple factor F, computes a suitable iterated pullback
  //                      of F (the factor p_B(F) in Gebhardt's paper).
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> MainPullback(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    list<CBraid::Factor<P>> ret;
    typename list<CBraid::Factor<P>>::iterator it = ret.end();

    CBraid::Braid<P> B2 = B;
    sint16 i;

    list<CBraid::Braid<P>> T = Trajectory(B);
    typename list<CBraid::Braid<P>>::reverse_iterator itb;

    CBraid::Factor<P> F2 = F;
    while (it == ret.end())
    {

      ret.push_back(F2);
      for (itb = T.rbegin(); itb != T.rend(); itb++)
        F2 = Pullback(*itb, F2);

      it = find(ret.begin(), ret.end(), F2);
    }

    typename list<CBraid::Factor<P>>::iterator it2 = it;

    sint16 l = 0;
    while (it2 != ret.end())
    {
      it2++;
      l++;
    }

    sint16 test;
    if (it == ret.begin())
    {
      test = 0;
    }
    else
      test = 1;

    it2 = ret.begin();
    while (test)
    {
      for (i = 0; i < l; i++)
      {
        it2++;
        if (it2 == it)
          test = 0;
      }
    }

    return *it2;
  }

  /////////////////////////////////////////////////////////////
  //
  //  MinUSS(B,F)  Given a braid B in its Ultra Summit Set (and in LCF),
  //               computes the minimal simple factor R such that
  //               F<R and B^R is in the Ultra Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> MinUSS(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Factor<P> F2 = MinSSS(B, F);

    list<CBraid::Factor<P>> ret = Returns(B, F2);
    typename list<CBraid::Factor<P>>::iterator it;

    for (it = ret.begin(); it != ret.end(); it++)
    {
      if (LeftMeet(F, *it) == F)
        return *it;
    }

    F2 = MainPullback(B, F);

    ret = Returns(B, F2);

    for (it = ret.begin(); it != ret.end(); it++)
    {
      if (LeftMeet(F, *it) == F)
        return *it;
    }
    cout << "Error in MinUSS.";
    exit(1);
  }

  /////////////////////////////////////////////////////////////
  //
  //  MinUSS(B)  Given a braid B in its Ultra Summit Set (and in LCF),
  //             computes the set of minimal simple factors R that
  //             B^R is in the Ultra Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Factor<P>> MinUSS(CBraid::Braid<P> B)
  {
    sint16 i, j, k, test;
    sint16 n = B.Index();
    sint16 *table = new sint16[n];

    list<CBraid::Factor<P>> Min;

    for (i = 0; i < n; i++)
      table[i] = 0;
    CBraid::Factor<P> F = CBraid::Factor<P>(n);
    for (i = 1; i < n; i++)
    {
      test = 1;
      F.Identity();
      k = F[i];
      F[i] = F[i + 1];
      F[i + 1] = k;

      F = MinUSS(B, F);

      for (j = 1; j < i; j++)
      {
        if (table[j - 1] && F[j] > F[j + 1])
          test = 0;
      }
      for (j = i + 1; j < n; j++)
      {
        if (F[j] > F[j + 1])
          test = 0;
      }
      if (test)
      {
        Min.push_back(F);
        table[i - 1] = 1;
      }
    }
    delete[] table;
    return Min;
  }

  /////////////////////////////////////////////////////////////
  //
  //  USS(B)  Given a braid B, computes its Ultra Summit Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<list<CBraid::Braid<P>>> USS(CBraid::Braid<P> B)
  {
    list<list<CBraid::Braid<P>>> uss;
    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index());
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itf, itf2;
    typename list<CBraid::Braid<P>>::iterator itb;

    CBraid::Braid<P> B2 = SendToUSS(B);
    list<CBraid::Braid<P>> T = Trajectory(B2);

    typename list<CBraid::Braid<P>>::reverse_iterator rit = T.rbegin();
    uss.push_back(Trajectory(Cycling(*rit)));

    B2 = ((!CBraid::Braid<P>(CBraid::Factor<P>(B.Index(), 1))) * (Cycling(*rit)) * CBraid::Factor<P>(B.Index(), 1)).MakeLCF();
    for (itb = (*uss.begin()).begin(); itb != (*uss.begin()).end(); itb++)
    {
      if (B2 == *itb)
        break;
    }

    if (itb == (*uss.begin()).end())
      uss.push_back(Trajectory(B2));

    typename list<list<CBraid::Braid<P>>>::iterator it = uss.begin(), it2;
    while (it != uss.end())
    {

      Min = MinUSS(*(*it).begin());

      for (itf = Min.begin(); itf != Min.end(); itf++)
      {
        F = *itf;

        B2 = ((!CBraid::Braid<P>(F)) * (*(*it).begin()) * F).MakeLCF();

        T = Trajectory(B2);
        for (itb = T.begin(); itb != T.end(); itb++)
        {
          for (it2 = uss.begin(); it2 != uss.end(); it2++)
          {
            if (*itb == *(*it2).begin())
              break;
          }
          if (it2 != uss.end())
            break;
        }

        if (itb == T.end())
        {
          uss.push_back(T);

          B2 = ((!CBraid::Braid<P>(CBraid::Factor<P>(B.Index(), 1))) * (*T.begin()) * CBraid::Factor<P>(B.Index(), 1)).MakeLCF();
          for (itb = T.begin(); itb != T.end(); itb++)
          {
            if (B2 == *itb)
              break;
          }

          if (itb == T.end())
            uss.push_back(Trajectory(B2));
        }
      }
      it++;
    }
    return uss;
  }

  /////////////////////////////////////////////////////////////
  //
  //  USS(B,mins,prev)  Given a braid B, computes its Ultra Summit Set,
  //                    and stores in the lists 'mins' and 'prev' the following data:
  //                    for each i, the first braid of the orbit i is obtained by
  //                    conjugation of the first element of the orbit prev[i]
  //                    by the simple element mins[i].
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<list<CBraid::Braid<P>>> USS(CBraid::Braid<P> B, list<CBraid::Factor<P>> &mins, list<sint16> &prev)
  {
    list<list<CBraid::Braid<P>>> uss;

    CBraid::Braid<P> B2 = SendToUSS(B);
    list<CBraid::Braid<P>> T = Trajectory(B2);
    typename list<CBraid::Braid<P>>::reverse_iterator rit = T.rbegin();
    uss.push_back(Trajectory(Cycling(*rit)));

    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index());
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itf, itf2;
    typename list<CBraid::Braid<P>>::iterator itb;

    sint16 current = 0;
    mins.clear();
    prev.clear();

    mins.push_back(CBraid::Factor<P>(B.Index(), 0));
    prev.push_back(1);

    typename list<list<CBraid::Braid<P>>>::iterator it = uss.begin(), it2;
    while (it != uss.end())
    {
      current++;

      Min = MinUSS(*(*it).begin());

      for (itf = Min.begin(); itf != Min.end(); itf++)
      {
        F = *itf;
        B2 = ((!CBraid::Braid<P>(F)) * (*(*it).begin()) * F).MakeLCF();
        T = Trajectory(B2);
        for (itb = T.begin(); itb != T.end(); itb++)
        {
          for (it2 = uss.begin(); it2 != uss.end(); it2++)
          {
            if (*itb == *(*it2).begin())
              break;
          }
          if (it2 != uss.end())
            break;
        }

        if (itb == T.end())
        {
          uss.push_back(T);
          mins.push_back(F);
          prev.push_back(current);
        }
      }
      it++;
    }

    return uss;
  }

  /////////////////////////////////////////////////////////////
  //
  //  TreePath(B,uss,mins,prev)  Computes a braid that conjugates
  //                             the first element in the Ultra Summit Set uss
  //                             to the braid B (which must be in the uss).
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> TreePath(CBraid::Braid<P> B, list<list<CBraid::Braid<P>>> &uss, list<CBraid::Factor<P>> &mins, list<sint16> &prev)
  {
    sint16 n = B.Index();
    CBraid::Braid<P> C = CBraid::Braid<P>(n);
    typename list<list<CBraid::Braid<P>>>::iterator it;
    typename list<CBraid::Braid<P>>::iterator itb, itb2;
    sint16 current = 0;
    list<sint16>::iterator itprev;
    typename list<CBraid::Factor<P>>::iterator itmins;
    sint16 i;

    if (CL(B) == 0)
      return CBraid::Braid<P>(n);

    for (it = uss.begin(); it != uss.end(); it++)
    {
      current++;
      for (itb = (*it).begin(); itb != (*it).end(); itb++)
      {
        if (*itb == B)
          break;
      }
      if (itb != (*it).end())
        break;
    }

    if (it == uss.end())
    {
      cout << "Error in TreePath" << endl;
      return 0;
    }

    for (itb2 = (*it).begin(); itb2 != itb; itb2++)
      C.RightMultiply((*(*itb2).FactorList.begin()).Flip(B.LeftDelta));

    while (current != 1)
    {
      itprev = prev.begin();
      itmins = mins.begin();
      for (i = 1; i < current; i++)
      {
        itprev++;
        itmins++;
      }
      C.LeftMultiply(*itmins);
      current = *itprev;
    }

    return C;
  }

  /////////////////////////////////////////////////////////////
  //
  //  AreConjugate(B1,B2,C)  Determines if the braids B1 and B2 are
  //                         conjugate, and computes a conjugating
  //                         element C.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  bool AreConjugate(CBraid::Braid<P> B1, CBraid::Braid<P> B2, CBraid::Braid<P> &C)
  {
    sint16 n = B1.Index();
    CBraid::Braid<P> C1 = CBraid::Braid<P>(n), C2 = CBraid::Braid<P>(n);

    CBraid::Braid<P> BT1 = SendToUSS(B1, C1), BT2 = SendToUSS(B2, C2);

    if (CL(BT1) != CL(BT2) || Sup(BT1) != Sup(BT2))
      return false;

    if (CL(BT1) == 0)
    {
      C = (C1 * (!C2)).MakeLCF();
      return true;
    }

    list<CBraid::Factor<P>> mins;
    list<sint16> prev;

    list<list<CBraid::Braid<P>>> uss = USS(BT1, mins, prev);

    typename list<list<CBraid::Braid<P>>>::iterator it;
    typename list<CBraid::Braid<P>>::iterator itb;
    sint16 current = 0;
    CBraid::Braid<P> D1 = CBraid::Braid<P>(n), D2 = CBraid::Braid<P>(n);

    for (it = uss.begin(); it != uss.end(); it++)
    {
      current++;
      D2 = CBraid::Braid<P>(n);
      for (itb = (*it).begin(); itb != (*it).end(); itb++)
      {
        if (*itb == BT2)
          break;
        D2 = D2 * ((*(*itb).FactorList.begin()).Flip((*itb).LeftDelta));
      }
      if (itb != (*it).end())
        break;
    }

    if (it == uss.end())
      return false;

    list<sint16>::iterator itprev;
    typename list<CBraid::Factor<P>>::iterator itmins;
    sint16 i;

    while (current != 1)
    {
      itprev = prev.begin();
      itmins = mins.begin();
      for (i = 1; i < current; i++)
      {
        itprev++;
        itmins++;
      }
      D1.LeftMultiply(*itmins);
      current = *itprev;
    }

    C = (C1 * D1 * D2 * (!C2)).MakeLCF();

    return true;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Centralizer(uss,mins,prev)  Computes the centralizer of the first
  //                              element in the Ultra Summit Set uss.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Braid<P>> Centralizer(list<list<CBraid::Braid<P>>> &uss, list<CBraid::Factor<P>> &mins, list<sint16> &prev)
  {
    CBraid::Braid<P> B = *(*uss.begin()).begin();
    sint16 n = B.Index();
    list<CBraid::Braid<P>> Cent;
    typename list<list<CBraid::Braid<P>>>::iterator it;
    typename list<CBraid::Braid<P>>::iterator itb;
    CBraid::Braid<P> C = CBraid::Braid<P>(n), D = CBraid::Braid<P>(n), E = CBraid::Braid<P>(n), B2 = CBraid::Braid<P>(n);
    sint16 cl = CL(B), sup = Sup(B), i;
    list<sint16> word;
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itMin;
    list<sint16>::iterator itprev;
    typename list<CBraid::Factor<P>>::iterator itmins;

    if (cl == 0 && sup % 2 == 0)
    {
      word.push_back(1);
      C = WordToBraid<P>(word, n);
      Cent.push_back(C);
      word.clear();

      for (i = 1; i < n; i++)
        word.push_back(i);
      C = WordToBraid<P>(word, n);
      Cent.push_back(C);

      return Cent;
    }

    if (cl == 0)
    {
      Min = MinUSS(B);
      for (itMin = Min.begin(); itMin != Min.end(); itMin++)
        Cent.push_back(CBraid::Braid<P>(*itMin));

      return Cent;
    }

    for (it = uss.begin(); it != uss.end(); it++)
    {
      D = TreePath(*(*it).begin(), uss, mins, prev);
      C = D;
      for (itb = (*it).begin(); itb != (*it).end(); itb++)
        C = C * ((*(*itb).FactorList.begin()).Flip(B.LeftDelta));
      C = C * (!D);
      C.MakeLCF();

      if (find(Cent.begin(), Cent.end(), C) == Cent.end())
        Cent.push_back(C);

      Min = MinUSS(*(*it).begin());
      for (itMin = Min.begin(); itMin != Min.end(); itMin++)
      {
        C = D * (*itMin);
        B2 = ((!CBraid::Braid<P>(*itMin)) * (*(*it).begin()) * (*itMin)).MakeLCF();
        E = TreePath(B2, uss, mins, prev);
        C = C * (!E);
        C.MakeLCF();

        if (C.CompareWithIdentity() == false && find(Cent.begin(), Cent.end(), C) == Cent.end())
          Cent.push_back(C);
      }
    }

    return Cent;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Centralizer(B)  Computes the centralizer of the braid B.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Braid<P>> Centralizer(CBraid::Braid<P> B)
  {
    sint16 n = B.Index();
    list<CBraid::Factor<P>> mins;
    list<sint16> prev;
    list<list<CBraid::Braid<P>>> uss = USS(B, mins, prev);

    list<CBraid::Braid<P>> Cent = Centralizer(uss, mins, prev);

    CBraid::Braid<P> C = CBraid::Braid<P>(n);
    SendToUSS(B, C);

    typename list<CBraid::Braid<P>>::iterator it;

    for (it = Cent.begin(); it != Cent.end(); it++)
    {
      (*it).LeftMultiply(C);
      (*it).RightMultiply(!C);
      (*it).MakeLCF();
    }

    return Cent;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Tableau(F,tab)  Computes the tableau associated to a
  //                  simple factor F, and stores it in tab.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  void Tableau(CBraid::Factor<P> F, sint16 **&tab)
  {
    sint16 i, j;
    sint16 n = F.Index();
    for (i = 0; i < n; i++)
    {
      tab[i][i] = F[i + 1];
    }
    for (j = 1; j <= n - 1; j++)
    {
      for (i = 0; i <= n - 1 - j; i++)
      {
        if (tab[i][i + j - 1] > tab[i + 1][i + j])
          tab[i][i + j] = tab[i][i + j - 1];
        else
          tab[i][i + j] = tab[i + 1][i + j];
      }
    }

    for (j = 1; j <= n - 1; j++)
    {
      for (i = j; i <= n - 1; i++)
      {
        if (tab[i - 1][i - j] < tab[i][i - j + 1])
          tab[i][i - j] = tab[i - 1][i - j];
        else
          tab[i][i - j] = tab[i][i - j + 1];
      }
    }
  }

  /////////////////////////////////////////////////////////////
  //
  //  Circles(B)  Determines if a braid B in LCF
  //              preserves a family of circles.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  bool Circles(CBraid::Braid<P> B)
  {
    sint16 j, k, t, d, n = B.Index();
    sint16 *disj = new sint16[n + 1];

    B.MakeLCF();

    sint16 cl = CL(B);
    sint16 delta, itype = 0;
    if (B.LeftDelta < 0)
      delta = -B.LeftDelta;
    else
      delta = B.LeftDelta;

    delta = delta % 2;

    sint16 ***tabarray = new sint16 **[cl + delta];
    typename list<CBraid::Factor<P>>::iterator it = B.FactorList.begin();

    for (j = 0; j < cl + delta; j++)
    {
      tabarray[j] = new sint16 *[n];
      for (k = 0; k < n; k++)
      {
        tabarray[j][k] = new sint16[n];
      }
      if (delta && j == 0)
        Tableau(CBraid::Factor<P>(n, 1), tabarray[j]);
      else
      {
        Tableau(*it, tabarray[j]);
        it++;
      }
    }

    sint16 *bkmove = new sint16[n];
    sint16 bk;
    for (j = 2; j < n; j++)
    {
      for (k = 1; k <= n - j + 1; k++)
      {
        bk = k;
        for (t = 0; t < cl + delta; t++)
        {
          if (tabarray[t][bk - 1][j + bk - 2] - tabarray[t][j + bk - 2][bk - 1] == j - 1)
            bk = tabarray[t][j + bk - 2][bk - 1];
          else
          {
            bk = 0;
            break;
          }
        }
        if (bk == k)
        {
          itype = 1;
          j = n + 1;
          break;
        }
        else if (bk - k < j && k - bk < j)
          bk = 0;

        bkmove[k] = bk;
      }
      for (k = 1; k <= n - j + 1; k++)
      {
        for (d = 1; d <= n; d++)
          disj[d] = 1;

        bk = k;
        while (bk)
        {
          if (bkmove[bk] == k)
          {
            itype = 1;
            k = n - j;
            j = n;
            break;
          }
          for (d = bk - j + 1; d <= bk + j - 1; d++)
          {
            if (d >= 1 && d <= n && d != k)
              disj[d] = 0;
          }
          bk = bkmove[bk];
          if (disj[bk] == 0)
            bk = 0;
        }
      }
    }

    if (itype)
      return true;
    else
      return false;
  }

  /////////////////////////////////////////////////////////////
  //
  //  ThurstonType(B)  Determines if a braid B is periodic (1),
  //                   reducible (2) or pseudo-Anosov (3).
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  int ThurstonType(CBraid::Braid<P> B)
  {
    sint16 i, n = B.Index();

    sint16 somereducible = 0, somePA = 0;

    B.MakeLCF();
    CBraid::Braid<P> pot = B;

    for (i = 0; i < n; i++)
    {
      if (CL(pot) == 0)
        return 1;
      pot = (pot * B).MakeLCF();
    }

    list<list<CBraid::Braid<P>>> uss = USS(B);
    typename list<list<CBraid::Braid<P>>>::iterator it;

    sint16 type = 3;

    for (it = uss.begin(); it != uss.end(); it++)
    {
      if (Circles(*(*it).begin()))
      {
        type = 2;
        somereducible = 1;
      }
      else
      {
        somePA = 1;
      }
    }

    if (somereducible && somePA)
      cout << "Not all elements in the USS preserve a family of circles!!!";

    return type;
  }

  /////////////////////////////////////////////////////////////
  //
  //  ThurstonType(uss)  Determines if the braids in the Ultra
  //                     Summit Set uss are periodic (1),
  //                     reducible (2) or pseudo-Anosov (3).
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  int ThurstonType(list<list<CBraid::Braid<P>>> &uss)
  {
    CBraid::Braid<P> B = *(*uss.begin()).begin();
    sint16 i, n = B.Index();

    sint16 somereducible = 0, somePA = 0;

    CBraid::Braid<P> pot = B;

    for (i = 0; i < n; i++)
    {
      if (CL(pot) == 0)
        return 1;
      pot = (pot * B).MakeLCF();
    }

    typename list<list<CBraid::Braid<P>>>::iterator it;

    sint16 type = 3;

    for (it = uss.begin(); it != uss.end(); it++)
    {
      if (Circles(*(*it).begin()))
      {
        type = 2;
        somereducible = 1;
      }
      else
      {
        somePA = 1;
      }
    }

    if (somereducible && somePA)
    {
      cout << "Not all elements in the USS of the braid " << endl;
      PrintBraidWord(*(*uss.begin()).begin());
      cout << endl
           << "preserve a family of circles!!!" << endl;
    }

    return type;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Rigidity(B)  Computes the rigidity of a braid B.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  sint16 Rigidity(CBraid::Braid<P> B)
  {
    CBraid::Braid<P> B2 = B.MakeLCF(), B3 = B2;
    sint16 cl = CL(B2), rigidity = 0;

    if (cl == 0)
      return rigidity;

    CBraid::Factor<P> F = (*B3.FactorList.begin()).Flip(B3.LeftDelta);
    B3 = B3 * F;
    B3.MakeLCF();

    typename list<CBraid::Factor<P>>::iterator it2, it3;

    it3 = B3.FactorList.begin();

    for (it2 = B2.FactorList.begin(); it2 != B2.FactorList.end(); it2++)
    {
      if (*it2 != *it3)
        break;
      rigidity++;
      it3++;
    }
    return rigidity;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Rigidity(uss)  Computes the maximal rigidity of a braid
  //                 in the Ultra Summit Set uss.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  sint16 Rigidity(list<list<CBraid::Braid<P>>> &uss)
  {
    typename list<list<CBraid::Braid<P>>>::iterator it;

    sint16 rigidity = 0, next, conjecture = 0;

    for (it = uss.begin(); it != uss.end(); it++)
    {
      if (it == uss.begin())
        rigidity = Rigidity(*(*it).begin());
      else
      {
        next = Rigidity(*(*it).begin());
        if (next != rigidity)
        {
          conjecture = 1;
          if (next > rigidity)
            rigidity = next;
        }
      }
    }
    if (conjecture)
    {
      cout << endl
           << "There are elements is the USS of" << endl;
      PrintBraidWord(*(*uss.begin()).begin());
      cout << endl
           << "with distinct rigidities!!!" << endl;
    }
    return rigidity;
  }

  /////////////////////////////////////////////////////////////
  //
  //  ReadIndex()   Asks to type the number of strands.
  //
  /////////////////////////////////////////////////////////////

  sint16 ReadIndex()
  {
    sint16 n;
    cout << endl
         << "Set the number of strands: ";
    std::cin >> n;
    std::cin.ignore();
    return n;
  }

  /////////////////////////////////////////////////////////////
  //
  //  ReadWord(n)   Asks to type a braid word on n strands,
  //               and returns the braid word.
  //
  /////////////////////////////////////////////////////////////

  list<sint16> ReadWord(sint16 n)
  {
    list<sint16> word;
    sint16 a;

    cout << endl
         << "Type a braid with " << n << " strands: "
         << "('" << n << "' = Delta)"
         << endl
         << endl;
    while (std::cin.peek() != '\n')
    {
      std::cin >> std::ws >> a;
      word.push_back(a);
    }
    std::cin.ignore();

    return word;
  }

  /////////////////////////////////////////////////////////////
  //
  //  ReadPower()   Asks to type the power to which the braid
  //                will be raised.
  //
  /////////////////////////////////////////////////////////////

  sint16 ReadPower()
  {
    sint16 power;
    cout << endl
         << "Raise it to power... ";
    std::cin >> power;
    std::cin.ignore();
    return power;
  }

  /////////////////////////////////////////////////////////////
  //
  //  RaisePower(B,k)   Raises the braid B to the power k.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> RaisePower(CBraid::Braid<P> B, sint16 k)
  {
    CBraid::Braid<P> original = B;
    sint16 i;
    if (k == 0)
      B.Identity();
    else if (k > 0)
    {
      for (i = 1; i < k; i++)
        B = B * original;
    }
    else
    {
      k = -k;
      original = !B;
      B = original;
      for (i = 1; i < k; i++)
        B = B * original;
    }

    return B;
  }

  /////////////////////////////////////////////////////////////
  //
  //  ReadFileName()   Asks to type the name of a file.
  //
  /////////////////////////////////////////////////////////////

  char *ReadFileName()
  {
    char *f = new char[30];
    cout << endl
         << "Type the name of the output file: ";
    std::cin.getline(f, 30);
    cout << endl;

    return f;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  //
  //  PrintUSS(uss,word,n,power,file,type,rigidity)   Prints the Ultra Summit Set
  //                                                  of the braid (word)^power to "file".
  //
  ////////////////////////////////////////////////////////////////////////////////////////

  template <class P>
  void PrintUSS(list<list<CBraid::Braid<P>>> &uss, list<sint16> word, sint16 n,
                sint16 power, char *file, sint16 type, sint16 rigidity)
  {
    std::ofstream f(file);

    sint16 orbits = 0;

    typename list<list<CBraid::Braid<P>>>::iterator it;
    typename list<CBraid::Braid<P>>::iterator oit, itb;
    CBraid::Braid<P> B2 = CBraid::Braid<P>(n);

    for (it = uss.begin(); it != uss.end(); it++)
      orbits++;

    sint16 *sizes = new sint16[orbits];
    sint16 size;
    sint16 i;

    it = uss.begin();
    for (i = 0; i < orbits; i++)
    {
      size = 0;
      for (oit = (*it).begin(); oit != (*it).end(); oit++)
      {
        size++;
      }
      sizes[i] = size;
      it++;
    }

    list<sint16>::iterator itw;

    f << "This file contains the Ultra Summit Set of the braid on " << n << " strands:" << endl
      << endl;
    if (power != 1)
      f << "( ";

    for (itw = word.begin(); itw != word.end(); itw++)
    {
      if (*itw == n)
        f << "D ";
      else if (*itw == -n)
        f << "-D ";
      else
        f << *itw << " ";
    }

    if (power != 1)
      f << ")^" << power;

    sint16 total;

    if (orbits == 1)
    {
      f << endl
        << endl
        << "It has 1 orbit, whose size is " << sizes[0] << "." << endl
        << endl;
      total = sizes[0];
    }
    else
    {
      f << endl
        << endl
        << "It has " << orbits << " orbits, whose sizes are: ";

      total = 0;

      for (i = 0; i < orbits; i++)
      {
        f << sizes[i];
        if (i == orbits - 1)
          f << "." << endl
            << endl;
        else
          f << ", ";
        total += sizes[i];
      }
    }

    f << "Total size:     " << total << endl
      << endl;

    if (type == 1)
      f << "Thurston type:  Periodic." << endl
        << endl;
    if (type == 2)
      f << "Thurston type:  Reducible." << endl
        << endl;
    if (type == 3)
      f << "Thurston type:  Pseudo-Anosov." << endl
        << endl;

    sint16 cl = CL(*(*uss.begin()).begin());

    f << "Rigidity:       ";

    if (rigidity == cl)
      f << "Rigid." << endl
        << endl;
    else if (rigidity == cl - 1 && cl != 1)
      f << "Subrigid." << endl
        << endl;
    else
      f << rigidity << "-rigid." << endl
        << endl;

    orbits = 1;

    for (it = uss.begin(); it != uss.end(); it++)
    {
      f << "-------------------" << endl
        << " Orbit number " << orbits++ << "        ";

      rigidity = Rigidity(*(*it).begin());
      if (rigidity == cl)
        f << "Rigid.";
      else if (rigidity == cl - 1 && cl != 1)
        f << "Subrigid.";
      else
        f << rigidity << "-rigid.";

      if (((!CBraid::Braid<P>(CBraid::Factor<P>(n, 1))) * (*(*it).begin()) * CBraid::Braid<P>(CBraid::Factor<P>(n, 1))).MakeLCF() == B2)
        f << "     Conjugate to the previous orbit by Delta.";
      else
      {
        for (itb = (*it).begin(); itb != (*it).end(); itb++)
        {
          if (((*(*it).begin()) * CBraid::Braid<P>(CBraid::Factor<P>(n, 1))).MakeLCF() == (CBraid::Braid<P>(CBraid::Factor<P>(n, 1)) * (*itb)).MakeLCF())
            break;
        }
        if (itb != (*it).end())
          f << "     Conjugate to itself by Delta.";
      }
      f << endl
        << "-------------------" << endl;

      B2 = *(*it).begin();

      size = 1;
      for (oit = (*it).begin(); oit != (*it).end(); oit++)
      {
        f << endl
          << std::setw(5) << size++;
        f << ":   ";
        f.close();
        PrintBraidWord(*oit, file);
        f.open(file, std::ios::app);
      }
      f << endl
        << endl
        << endl;
    }
  }

  /////////////////////////////////////////////////////////////
  //
  //   FileName(iteration,max_iteration,type,orbit,rigidity,cl)  Creates the file name
  //                                                       corresponding to the given data.
  //
  /////////////////////////////////////////////////////////////

  char *FileName(sint16 iteration, sint16 max_iteration, sint16 type, sint16 orbit, sint16 rigidity, sint16 cl)
  {
    char *file = new char[30];
    sint16 i, j;

    if (type == 1)
    {
      file[0] = 'p';
      file[1] = 'e';
    }
    else if (type == 2)
    {
      file[0] = 'r';
      file[1] = 'e';
    }
    else
    {
      file[0] = 'p';
      file[1] = 'a';
    }

    file[2] = '_';

    if (rigidity == cl)
      file[3] = 'R';
    else if (rigidity == cl - 1 && cl != 1)
      file[3] = 'S';
    else
      file[3] = '0' + rigidity;

    file[4] = '_';

    if (orbit > 9)
      file[5] = 'M';
    else
      file[5] = '0' + orbit;

    file[6] = '_';

    sint16 digits = 1;
    j = 10;
    while (max_iteration / j >= 1)
    {
      j = j * 10;
      digits++;
    }

    i = 1;
    for (j = 1; j < digits; j++)
      i = i * 10;

    sint16 partial_iteration = iteration;
    for (j = 7; j <= digits + 6; j++)
    {
      file[j] = '0' + (partial_iteration / i);
      partial_iteration -= (partial_iteration / i) * i;
      i = i / 10;
    }

    file[digits + 7] = '.';
    file[digits + 8] = 't';
    file[digits + 9] = 'x';
    file[digits + 10] = 't';
    file[digits + 11] = 0;

    return file;
  }

  //////////////////---------------

  ///////////////////////////////////////////////////////
  //
  //  Reverse(B)  computes the revese of a braid B,
  //              that is, B written backwards.
  //              B must be given in left canonical form.
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> Reverse(CBraid::Braid<P> B)
  {
    sint16 i, l = CL(B);
    CBraid::Braid<P> B2 = CBraid::Braid<P>(B.Index());
    B2.RightDelta = B.LeftDelta;

    for (i = 0; i < l; i++)
    {
      B2.FactorList.push_front(!B.FactorList.front());
      B.FactorList.pop_front();
    }
    B2.MakeLCF();
    return B2;
  }

  /////////////////////////////////////////////////////////////
  //
  //  RightMeet(B1,B2)  Given two braids B1 and B2, computes
  //                    their right gcd. That is, the greatest braid
  //                    B such that B1>B and B2>B.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> RightMeet(CBraid::Braid<P> B1, CBraid::Braid<P> B2)
  {
    return Reverse(LeftMeet(Reverse(B1), Reverse(B2)));
  }

  /////////////////////////////////////////////////////////////
  //
  //  LeftJoin(B1,B2)  Given two braids B1 and B2, computes
  //                    their left lcm. That is, the smallest braid
  //                    B such that B1<B and B2<B.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> LeftJoin(CBraid::Braid<P> B1, CBraid::Braid<P> B2)
  {
    sint16 n = B1.Index(), shift = 0;
    CBraid::Braid<P> B = CBraid::Braid<P>(n);
    CBraid::Factor<P> F2 = CBraid::Factor<P>(n, 0), F = CBraid::Factor<P>(n, 1);

    B1.MakeLCF();
    B2.MakeLCF();

    shift = (B1.LeftDelta < B2.LeftDelta ? B1.LeftDelta : B2.LeftDelta);

    B1.LeftDelta -= shift;
    B2.LeftDelta -= shift;

    B = B1;

    while (!B2.CompareWithIdentity())
    {
      if (B2.LeftDelta > 0)
        F2 = CBraid::Factor<P>(n, 1);
      else if (CL(B2) == 0)
        F2 = CBraid::Factor<P>(n, 0);
      else
        F2 = B2.FactorList.front();

      F = Remainder(B1, F2);

      B.RightMultiply(F);
      B1.RightMultiply(F);

      ////////////// Multiply B1 from the left by F2^{-1}.

      B1.LeftDelta--;
      B1.FactorList.push_front((~F2).Flip(B1.LeftDelta));
      B1.MakeLCF();

      ////////////// Multiply B2 from the left by F2^{-1}.

      B2.LeftDelta--;
      B2.FactorList.push_front((~F2).Flip(B2.LeftDelta));
      B2.MakeLCF();
    }

    B.MakeLCF();
    B.LeftDelta += shift;
    return B;
  }
  /////////////////////////////////////////////////////////////
  //
  //  RightJoin(B1,B2)  Given two braids B1 and B2, computes
  //                    their right lcm. That is, the smallest braid
  //                    B such that B>B1 and B>B2.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> RightJoin(CBraid::Braid<P> B1, CBraid::Braid<P> B2)
  {
    return Reverse(LeftJoin(Reverse(B1), Reverse(B2)));
  }

  ///////////////////////////////////////////////////////
  //
  //  InitialFactor(B)  computes the initial factor of a braid B,
  //                    given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> InitialFactor(CBraid::Braid<P> B)
  {
    sint16 n = B.Index();
    CBraid::Factor<P> F = CBraid::Factor<P>(n, 0);

    if (CL(B) > 0)
      F = (B.FactorList.front()).Flip(-B.LeftDelta);

    return F;
  }

  ///////////////////////////////////////////////////////
  //
  //  PreferredPrefix(B)  computes the preferred prefix of a braid B,
  //                      given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> PreferredPrefix(CBraid::Braid<P> B)
  {
    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index(), 0);

    if (CL(B) > 0)
      F = LeftMeet(InitialFactor(B), ~B.FactorList.back());

    return F;
  }

  ///////////////////////////////////////////////////////
  //
  //  Sliding(B)  computes the cyclic sliding of a braid B,
  //                given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> Sliding(CBraid::Braid<P> B)
  {
    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index());

    if (CL(B) == 0)
      return B;

    F = PreferredPrefix(B);
    B.FactorList.front() = (!(F.Flip(B.LeftDelta))) * B.FactorList.front();
    B.FactorList.push_back(F);

    return B.MakeLCF();
  }

  ///////////////////////////////////////////////////////
  //
  //  PreferredSuffix(B)  computes the preferred suffix of a braid B,
  //                      given in Left Canonical Form
  //
  ///////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> PreferredSuffix(CBraid::Braid<P> B)
  {
    return !(PreferredPrefix(Reverse(B)));
  }

  /////////////////////////////////////////////////////////////
  //
  //  Trajectory_Sliding(B)  Computes the trajectory under cyclic sliding
  //                         of a braid B, that is, a list containing eta^k(B),
  //                         for k=0,1,... until the first repetition.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Braid<P>> Trajectory_Sliding(CBraid::Braid<P> B)
  {
    list<CBraid::Braid<P>> p;

    while (find(p.begin(), p.end(), B) == p.end())
    {
      p.push_back(B);
      B = Sliding(B);
    }
    return p;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Trajectory_Sliding(B,C,d)  Computes the trajectory of a braid B for cyclic sliding,
  //                     a braid C that conjugates B to the
  //                     first element of a closed orbit under sliding,
  //                     and the number d of slidings needed to reach that element
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Braid<P>> Trajectory_Sliding(CBraid::Braid<P> B, CBraid::Braid<P> &C, sint16 &d)
  {
    list<CBraid::Braid<P>> p;
    typename list<CBraid::Braid<P>>::iterator it;

    sint16 n = B.Index();
    C = CBraid::Braid<P>(n);
    d = 0;

    while (find(p.begin(), p.end(), B) == p.end())
    {
      p.push_back(B);
      C.RightMultiply(PreferredPrefix(B));
      B = Sliding(B);
      d++;
    }

    CBraid::Braid<P> B2 = CBraid::Braid<P>(n);
    CBraid::Braid<P> C2 = CBraid::Braid<P>(n);

    C2.RightMultiply(PreferredPrefix(B));
    B2 = Sliding(B);
    d--;
    while (B2 != B)
    {
      C2.RightMultiply(PreferredPrefix(B2));
      B2 = Sliding(B2);
      d--;
    }
    C.RightMultiply(!C2);
    C.MakeLCF();

    return p;
  }

  /////////////////////////////////////////////////////////////
  //
  //  SendToSC(B)  Computes a braid conjugate to B that
  //                belongs to its Sliding Circuits Set.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> SendToSC(CBraid::Braid<P> B)
  {
    list<CBraid::Braid<P>> T = Trajectory_Sliding(B);
    return Sliding(T.back());
  }

  /////////////////////////////////////////////////////////////
  //
  //  SendToSC(B,C)  Computes a braid conjugate to B that
  //                  belongs to its Sliding Circuits Set, and a braid
  //                  C that conjugates B to the result.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Braid<P> SendToSC(CBraid::Braid<P> B, CBraid::Braid<P> &C)
  {
    sint16 d;
    list<CBraid::Braid<P>> T = Trajectory_Sliding(B, C, d);

    return Sliding(T.back());
  }

  /////////////////////////////////////////////////////////////
  //
  //  Transport_Sliding(B,F)   Given a braid B (in its SC and in LCF), and a simple factor
  //                   F such that B^F is in its SSS, computes the transport of F for sliding.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> Transport_Sliding(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Braid<P> B2 = ((!CBraid::Braid<P>(F)) * B * F).MakeLCF();
    CBraid::Braid<P> B3 = ((!CBraid::Braid<P>(PreferredPrefix(B))) * F * (PreferredPrefix(B2))).MakeLCF();

    CBraid::Factor<P> F2 = CBraid::Factor<P>(B2.Index(), 0);

    if (CL(B3) > 0)
      F2 = B3.FactorList.front();
    else if (B3.LeftDelta == 1)
      F2 = CBraid::Factor<P>(B2.Index(), 1);

    return F2;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Returns_Sliding(B,F)   Given a braid B (in its SC and in LCF), and a simple factor
  //                 F such that B^F is in its SSS, computes the iterated transports
  //                 of F for sliding that send B to an element in the circuit of B^F.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Factor<P>> Returns_Sliding(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    list<CBraid::Factor<P>> ret;
    typename list<CBraid::Factor<P>>::iterator it = ret.end();
    CBraid::Braid<P> B1 = B;
    sint16 i, N = 1;

    B1 = Sliding(B1);
    while (B1 != B)
    {
      N++;
      B1 = Sliding(B1);
    }

    while (it == ret.end())
    {
      ret.push_back(F);

      B1 = B;
      for (i = 0; i < N; i++)
      {
        F = Transport_Sliding(B1, F);
        B1 = Sliding(B1);
      }

      it = find(ret.begin(), ret.end(), F);
    }

    while (it != ret.begin())
      ret.pop_front();

    return ret;
  }

  /////////////////////////////////////////////////////////////
  //
  //  Pullback_Sliding(B,F)   Given a braid B (in its SC and in LCF), and a
  //                          simple factor F such that B^F is super summit,
  //                          computes the pullback of F at s(B) for sliding.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> Pullback_Sliding(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {

    CBraid::Braid<P> B2 = CBraid::Braid<P>(B.Index());
    B2 = (CBraid::Braid<P>(PreferredPrefix(B)) * F).MakeLCF();

    CBraid::Braid<P> B3 = CBraid::Braid<P>(B.Index());
    B3 = (!CBraid::Braid<P>(F) * (Sliding(B) * F)).MakeLCF();

    CBraid::Factor<P> F2 = CBraid::Factor<P>(B.Index());
    F2 = PreferredSuffix(B3);

    CBraid::Braid<P> C = CBraid::Braid<P>(B.Index());
    C = RightMeet(B2, CBraid::Braid<P>(F2));

    B2 = (B2 * (!C)).MakeLCF();

    if (B2.CompareWithIdentity())
      return CBraid::Factor<P>(B.Index(), 0);
    else if (CL(B2) == 0)
      return CBraid::Factor<P>(B.Index(), 1);
    else
      return B2.FactorList.front();
  }

  /////////////////////////////////////////////////////////////
  //
  //  MainPullback_Sliding(B,F)   Given a braid B (in its SC and in LCF), and a
  //                      simple factor F, computes the first repeated iterated
  //                      pullback for cyclic sliding of F.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> MainPullback_Sliding(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {

    list<CBraid::Factor<P>> ret;
    typename list<CBraid::Factor<P>>::iterator it = ret.end();

    CBraid::Braid<P> B2 = B;

    list<CBraid::Braid<P>> T = Trajectory_Sliding(B);
    typename list<CBraid::Braid<P>>::reverse_iterator itb;

    if (F.CompareWithDelta())
      return F;

    CBraid::Factor<P> F2 = F;
    while (it == ret.end())
    {

      ret.push_back(F2);
      for (itb = T.rbegin(); itb != T.rend(); itb++)
        F2 = Pullback_Sliding(*itb, F2);

      it = find(ret.begin(), ret.end(), F2);
    }

    return *it;
  }

  // Mar�a Cumplido Cabello

  /////////////////////////////////////////////////////////////
  //
  //  MinSC(B,F)  Given a braid B in its Set of Sliding Circuits (and in LCF),
  //               computes the minimal simple factor R such that
  //               F<R and B^R is in the Set of Sliding Circuits.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  CBraid::Factor<P> MinSC(CBraid::Braid<P> B, CBraid::Factor<P> F)
  {
    CBraid::Factor<P> F2 = MinSSS(B, F);

    list<CBraid::Factor<P>> ret = Returns_Sliding(B, F2);
    typename list<CBraid::Factor<P>>::iterator it;

    for (it = ret.begin(); it != ret.end(); it++)
    {
      if (LeftMeet(F, *it) == F)
        return *it;
    }

    F2 = MainPullback_Sliding(B, F);

    ret = Returns_Sliding(B, F2);

    for (it = ret.begin(); it != ret.end(); it++)
    {
      if (LeftMeet(F, *it) == F)
        return *it;
    }

    return CBraid::Factor<P>(B.Index(), 1);
  }

  /////////////////////////////////////////////////////////////
  //
  //  MinSC(B)  Given a braid B in its Set of Sliding Circuits (and in LCF),
  //             computes the set of minimal simple factors R that
  //             B^R is in the Set of Sliding Circuits.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<CBraid::Factor<P>> MinSC(CBraid::Braid<P> B)
  {
    sint16 i, j, k, test;
    sint16 n = B.Index();
    sint16 *table = new sint16[n];

    list<CBraid::Factor<P>> Min;

    for (i = 0; i < n; i++)
      table[i] = 0;
    CBraid::Factor<P> F = CBraid::Factor<P>(n);
    for (i = 1; i < n; i++)
    {
      test = 1;
      F.Identity();
      k = F[i];
      F[i] = F[i + 1];
      F[i + 1] = k;

      F = MinSC(B, F);

      for (j = 1; j < i; j++)
      {
        if (table[j - 1] && F[j] > F[j + 1])
          test = 0;
      }
      for (j = i + 1; j < n; j++)
      {
        if (F[j] > F[j + 1])
          test = 0;
      }
      if (test)
      {
        Min.push_back(F);
        table[i - 1] = 1;
      }
    }
    delete[] table;
    return Min;
  }

  /////////////////////////////////////////////////////////////
  //
  //  SC(B)  Given a braid B, computes its Set of Cyclic Slidings.
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<list<CBraid::Braid<P>>> SC(CBraid::Braid<P> B)
  {
    list<list<CBraid::Braid<P>>> sc;
    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index());
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itf, itf2;
    typename list<CBraid::Braid<P>>::iterator itb;

    CBraid::Braid<P> B2 = SendToSC(B);
    list<CBraid::Braid<P>> T = Trajectory_Sliding(B2);

    sc.push_back(Trajectory_Sliding(B2));

    B2 = ((!CBraid::Braid<P>(CBraid::Factor<P>(B.Index(), 1))) * (B2)*CBraid::Factor<P>(B.Index(), 1)).MakeLCF();

    for (itb = (*sc.begin()).begin(); itb != (*sc.begin()).end(); itb++)
    {
      if (B2 == *itb)
        break;
    }

    if (itb == (*sc.begin()).end())
      sc.push_back(Trajectory_Sliding(B2));

    typename list<list<CBraid::Braid<P>>>::iterator it = sc.begin(), it2;

    while (it != sc.end())
    {

      Min = MinSC(*(*it).begin());

      for (itf = Min.begin(); itf != Min.end(); itf++)
      {
        F = *itf;
        B2 = ((!CBraid::Braid<P>(F)) * (*(*it).begin()) * F).MakeLCF();
        T = Trajectory_Sliding(B2);
        for (itb = T.begin(); itb != T.end(); itb++)
        {
          for (it2 = sc.begin(); it2 != sc.end(); it2++)
          {
            if (*itb == *(*it2).begin())
              break;
          }
          if (it2 != sc.end())
            break;
        }

        if (itb == T.end())
        {
          sc.push_back(T);
          B2 = ((!CBraid::Braid<P>(CBraid::Factor<P>(B.Index(), 1))) * (*T.begin()) * CBraid::Factor<P>(B.Index(), 1)).MakeLCF();
          for (itb = T.begin(); itb != T.end(); itb++)
          {
            if (B2 == *itb)
              break;
          }
          if (itb == T.end())
            sc.push_back(Trajectory_Sliding(B2));
        }
      }
      it++;
    }
    return sc;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  //
  //  PrintSC(sc,word,n,power,file,type)   Prints the Set of Sliding Circuits
  //                                       of the braid (word)^power to "file".
  //
  ////////////////////////////////////////////////////////////////////////////////////////

  template <class P>
  void PrintSC(list<list<CBraid::Braid<P>>> &sc, list<sint16> word, sint16 n,
               sint16 power, char *file, sint16 type)
  {
    std::ofstream f(file);

    sint16 orbits = 0;

    typename list<list<CBraid::Braid<P>>>::iterator it;
    typename list<CBraid::Braid<P>>::iterator oit, itb;
    CBraid::Braid<P> B2 = CBraid::Braid<P>(n);

    for (it = sc.begin(); it != sc.end(); it++)
      orbits++;

    sint16 *sizes = new sint16[orbits];
    sint16 size;
    sint16 i;

    it = sc.begin();
    for (i = 0; i < orbits; i++)
    {
      size = 0;
      for (oit = (*it).begin(); oit != (*it).end(); oit++)
      {
        size++;
      }
      sizes[i] = size;
      it++;
    }

    list<sint16>::iterator itw;

    f << "This file contains the Set of Sliding Circuits of the braid on " << n << " strands:" << endl
      << endl;
    if (power != 1)
      f << "( ";

    for (itw = word.begin(); itw != word.end(); itw++)
    {
      if (*itw == n)
        f << "D ";
      else if (*itw == -n)
        f << "-D ";
      else
        f << *itw << " ";
    }

    if (power != 1)
      f << ")^" << power;

    sint16 total;

    if (orbits == 1)
    {
      f << endl
        << endl
        << "It has 1 circuit, whose size is " << sizes[0] << "." << endl
        << endl;
      total = sizes[0];
    }
    else
    {
      f << endl
        << endl
        << "It has " << orbits << " circuits, whose sizes are: ";

      total = 0;

      for (i = 0; i < orbits; i++)
      {
        f << sizes[i];
        if (i == orbits - 1)
          f << "." << endl
            << endl;
        else
          f << ", ";
        total += sizes[i];
      }
    }

    f << "Total size:     " << total << endl
      << endl;

    if (type == 1)
      f << "Thurston type:  Periodic." << endl
        << endl;
    if (type == 2)
      f << "Thurston type:  Reducible." << endl
        << endl;
    if (type == 3)
      f << "Thurston type:  Pseudo-Anosov." << endl
        << endl;

    orbits = 1;

    for (it = sc.begin(); it != sc.end(); it++)
    {
      f << "-------------------" << endl
        << " Circuit number " << orbits++ << "        ";

      if (((!CBraid::Braid<P>(CBraid::Factor<P>(n, 1))) * (*(*it).begin()) * CBraid::Braid<P>(CBraid::Factor<P>(n, 1))).MakeLCF() == B2)
        f << "     Conjugate to the previous circuit by Delta.";
      else
      {
        for (itb = (*it).begin(); itb != (*it).end(); itb++)
        {
          if (((*(*it).begin()) * CBraid::Braid<P>(CBraid::Factor<P>(n, 1))).MakeLCF() == (CBraid::Braid<P>(CBraid::Factor<P>(n, 1)) * (*itb)).MakeLCF())
            break;
        }
        if (itb != (*it).end())
          f << "     Conjugate to itself by Delta.";
      }
      f << endl
        << "-------------------" << endl;

      B2 = *(*it).begin();

      size = 1;
      for (oit = (*it).begin(); oit != (*it).end(); oit++)
      {
        f << endl
          << std::setw(5) << size++;
        f << ":   ";
        f.close();
        PrintBraidWord(*oit, file);
        f.open(file, std::ios::app);
      }
      f << endl
        << endl
        << endl;
    }
  }

  /////////////////////////////////////////////////////////////
  //
  //  SC(B,mins,prev)  Given a braid B, computes its Set of Sliding Circuits,
  //                    and stores in the lists 'mins' and 'prev' the following data:
  //                    for each i, the first braid of the orbit i is obtained by
  //                    conjugation of the first element of the orbit prev[i]
  //                    by the simple element mins[i].
  //
  /////////////////////////////////////////////////////////////

  template <class P>
  list<list<CBraid::Braid<P>>> SC(CBraid::Braid<P> B, list<CBraid::Factor<P>> &mins, list<sint16> &prev)
  {
    list<list<CBraid::Braid<P>>> sc;

    CBraid::Braid<P> B2 = SendToSC(B);
    list<CBraid::Braid<P>> T = Trajectory_Sliding(B2);
    sc.push_back(Trajectory_Sliding(B2));

    CBraid::Factor<P> F = CBraid::Factor<P>(B.Index());
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itf, itf2;
    typename list<CBraid::Braid<P>>::iterator itb;

    sint16 current = 0;
    mins.clear();
    prev.clear();

    mins.push_back(CBraid::Factor<P>(B.Index(), 0));
    prev.push_back(1);

    typename list<list<CBraid::Braid<P>>>::iterator it = sc.begin(), it2;
    while (it != sc.end())
    {
      current++;

      Min = MinSC(*(*it).begin());

      for (itf = Min.begin(); itf != Min.end(); itf++)
      {
        F = *itf;
        B2 = ((!CBraid::Braid<P>(F)) * (*(*it).begin()) * F).MakeLCF();
        T = Trajectory_Sliding(B2);
        for (itb = T.begin(); itb != T.end(); itb++)
        {
          for (it2 = sc.begin(); it2 != sc.end(); it2++)
          {
            if (*itb == *(*it2).begin())
              break;
          }
          if (it2 != sc.end())
            break;
        }

        if (itb == T.end())
        {
          sc.push_back(T);
          mins.push_back(F);
          prev.push_back(current);
        }
      }
      it++;
    }

    return sc;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  //  AreConjugateSC(B1,B2,C)  Determines if the braids B1 and B2 are
  //                           conjugate by testing their set of sliding circuits,
  //                           and computes a conjugating element C.
  //
  //////////////////////////////////////////////////////////////////////////////////

  template <class P>
  bool AreConjugateSC(CBraid::Braid<P> B1, CBraid::Braid<P> B2, CBraid::Braid<P> &C)
  {
    sint16 n = B1.Index();
    CBraid::Braid<P> C1 = CBraid::Braid<P>(n), C2 = CBraid::Braid<P>(n);

    CBraid::Braid<P> BT1 = SendToSC(B1, C1), BT2 = SendToSC(B2, C2);

    if (CL(BT1) != CL(BT2) || Sup(BT1) != Sup(BT2))
      return false;

    if (CL(BT1) == 0)
    {
      C = (C1 * (!C2)).MakeLCF();
      return true;
    }

    list<CBraid::Factor<P>> mins;
    list<sint16> prev;

    list<list<CBraid::Braid<P>>> sc = SC(BT1, mins, prev);

    typename list<list<CBraid::Braid<P>>>::iterator it;
    typename list<CBraid::Braid<P>>::iterator itb;
    sint16 current = 0;
    CBraid::Braid<P> D1 = CBraid::Braid<P>(n), D2 = CBraid::Braid<P>(n);

    for (it = sc.begin(); it != sc.end(); it++)
    {
      current++;
      D2 = CBraid::Braid<P>(n);
      for (itb = (*it).begin(); itb != (*it).end(); itb++)
      {
        if (*itb == BT2)
          break;
        D2 = D2 * PreferredPrefix(*itb);
      }
      if (itb != (*it).end())
        break;
    }

    if (it == sc.end())
      return false;

    list<sint16>::iterator itprev;
    typename list<CBraid::Factor<P>>::iterator itmins;
    sint16 i;

    while (current != 1)
    {
      itprev = prev.begin();
      itmins = mins.begin();
      for (i = 1; i < current; i++)
      {
        itprev++;
        itmins++;
      }
      D1.LeftMultiply(*itmins);
      current = *itprev;
    }

    C = (C1 * D1 * D2 * (!C2)).MakeLCF();

    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  //  AreConjugateSC2(B1,B2,C)  Determines if the braids B1 and B2 are
  //                           conjugate by the procedure of contruct SC(B1),
  //                           and computes a conjugating element C.
  //
  //////////////////////////////////////////////////////////////////////////////////

  template <class P>
  bool AreConjugateSC2(CBraid::Braid<P> B1, CBraid::Braid<P> B2, CBraid::Braid<P> &C)
  {
    list<list<CBraid::Braid<P>>> sc;
    sint16 n = B1.Index();
    CBraid::Factor<P> F = CBraid::Factor<P>(n);
    list<CBraid::Factor<P>> Min;
    typename list<CBraid::Factor<P>>::iterator itf, itf2;
    typename list<CBraid::Braid<P>>::iterator itb;

    CBraid::Braid<P> C1 = CBraid::Braid<P>(n), C2 = CBraid::Braid<P>(n);
    CBraid::Braid<P> BT1 = SendToSC(B1, C1), BT2 = SendToSC(B2, C2);

    CBraid::Braid<P> D1 = CBraid::Braid<P>(n), D2 = CBraid::Braid<P>(n);

    list<CBraid::Factor<P>> mins;
    list<sint16> prev;

    sint16 current = 0;
    sint16 current2 = 0;
    mins.clear();
    prev.clear();

    list<sint16>::iterator itprev;
    typename list<CBraid::Factor<P>>::iterator itmins;
    sint16 i;

    mins.push_back(CBraid::Factor<P>(B1.Index(), 0));
    prev.push_back(1);

    list<CBraid::Braid<P>> T = Trajectory_Sliding(BT1);

    sc.push_back(Trajectory_Sliding(BT1));

    typename list<list<CBraid::Braid<P>>>::iterator it = sc.begin(), it2;

    while (it != sc.end())
    {
      current++;
      current2++;
      D2 = CBraid::Braid<P>(n);
      for (itb = (*it).begin(); itb != (*it).end(); itb++)
      {
        if (*itb == BT2)
        {
          while (current2 != 1)
          {
            itprev = prev.begin();
            itmins = mins.begin();
            for (i = 1; i < current2; i++)
            {
              itprev++;
              itmins++;
            }
            D1.LeftMultiply(*itmins);
            current2 = *itprev;
          }
          C = (C1 * D1 * D2 * (!C2)).MakeLCF();
          return true;
        }
        D2 = D2 * PreferredPrefix(*itb);
      }

      Min = MinSC(*(*it).begin());

      for (itf = Min.begin(); itf != Min.end(); itf++)
      {
        F = *itf;
        BT1 = ((!CBraid::Braid<P>(F)) * (*(*it).begin()) * F).MakeLCF();

        T = Trajectory_Sliding(BT1);
        for (itb = T.begin(); itb != T.end(); itb++)
        {
          for (it2 = sc.begin(); it2 != sc.end(); it2++)
          {
            if (*itb == *(*it2).begin())
              break;
          }
          if (it2 != sc.end())
            break;
        }

        if (itb == T.end())
        {
          sc.push_back(T);
          mins.push_back(F);
          prev.push_back(current);
        }
      }
      it++;
    }
    return false;
  }

}

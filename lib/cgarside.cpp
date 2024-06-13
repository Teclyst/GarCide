#include "cgarside.h"
#include <iostream>
#include <ostream>

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

  void foo()
  {

    ArtinBraidUnderlying bar = ArtinBraidUnderlying(3);

    bar.Delta();

    ArtinBraidFactor id = ArtinBraidFactor(bar.Copy());

    id.Neutral();

    id.Randomize();

    id.Print(std::cout);
  }

}
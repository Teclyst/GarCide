#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

#include "artin_braid.h"
#include "band_braid.h"
#include "octaedral_braid.h"
#include "dihedral_braid.h"
#include "garsidesc.h"
#include "garsideuss.h"
#include "cbraid.h"
#include "braiding.h"

#include "optarg.h"
#include "timecounter.h"

int Index = 24;
int CLength = 1;
int Count = 1;
int RandomSeed;

CBraid::ArtinFactor ToCBraidFactor(const CGarside::ArtinBraidFactor &f) {
  CBraid::ArtinFactor cf(f.GetParameter());
  CGarside::ArtinBraidUnderlying u = f.GetUnderlying();
  for (int i = 1; i <= f.GetParameter(); i++) {
    cf.At(i) = u.At(i);
  }
  return cf;
}

CBraid::ArtinBraid ToCBraid(const CGarside::ArtinBraid &b) {
  CBraid::ArtinBraid cb(b.GetParameter());
  cb.LeftDelta = b.Inf();
  for (CGarside::ArtinBraid::ConstFactorItr it = b.FactorList.begin(); it != b.FactorList.end(); it++) {
    cb.FactorList.push_back(ToCBraidFactor(*it));
  }
  return cb;
}

int main()
{ 
  srand(RandomSeed);

  CGarside::ArtinBraid b(Index);

  b.Randomize(CLength);
  b.Normalize();

  CBraid::ArtinBraid cb = ToCBraid(b);

  TimeCounter t;

  t.Start();

  std::list<std::list<CBraid::ArtinBraid>> cscs = Braiding::SC(cb);

  t.Stop();

  std::cout << "Computing (Braiding) the SCS took " << t.IntervalSec() * 1000 << " ms." << std::endl;
  std::cout << "It has " << cscs.size() << " orbits, and " << cscs.front().size() * cscs.size()  << " elements." << std::endl;

  //t.Start();

  //SSS::SuperSummitSet<CGarside::ArtinBraid> sss = SSS::SSS(b);

  //t.Stop();

  //std::cout << "Computing the SSS, of size " << sss.Card() << ", took " << t.IntervalSec() * 1000 << " ms." << std::endl;

  t.Start();

  SC::SlidingCircuitSet<CGarside::ArtinBraid> scs = SC::SCS(b);

  t.Stop();

  //scs.Print(std::cout);

  std::cout << "Computing the SCS took " << t.IntervalSec() * 1000 << " ms." << std::endl;
  std::cout << "It has " << scs.NumberOfOrbits() << " orbits, and " << scs.Card()  << " elements." << std::endl;

  return 0; 
}
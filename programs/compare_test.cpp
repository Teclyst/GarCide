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

int Index = 20;
int CLength = 256;
int Count = 1;
int RandomSeed = 8577;

int main()
{ 
  srand(RandomSeed);

  CGarside::ArtinBraid b(Index);
  CBraid::ArtinBraid cb(Index);
  cb.Randomize(CLength);
  cb.MakeLCF();

  b.Randomize(CLength);
  b.Normalize();

  //b.Print(std::cout);
  //std::cout << std::endl;

  TimeCounter t;

  //t.Start();
//
  //Braiding::SC(cb);
//
  //t.Stop();
//
  //std::cout << "Computing (Braiding) the SCS took " << t.IntervalSec() * 1000 << " ms." << std::endl;

  //t.Start();

  //SSS::SuperSummitSet<CGarside::ArtinBraid> sss = SSS::SSS(b);

  //t.Stop();

  //std::cout << "Computing the SSS, of size " << sss.Card() << ", took " << t.IntervalSec() * 1000 << " ms." << std::endl;

  t.Start();

  SC::SlidingCircuitSet<CGarside::ArtinBraid> scs = SC::SCS(b);

  t.Stop();

  //scs.Print(std::cout);

  std::cout << "Computing the SCS took " << t.IntervalSec() * 1000 << " ms." << std::endl;

  return 0; 
}
#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

#include "artin_braid.h"
#include "band_braid.h"
#include "octaedral_braid.h"
#include "dihedral_braid.h"
#include "dual_complex_reflection.h"
#include "garsidesc.h"
#include "garsideuss.h"
#include "cbraid.h"
#include "braiding.h"
#include <omp.h>
#include "optarg.h"
#include "timecounter.h"

int Index = 4;
int CLength = 256;
int Count = 1;
int RandomSeed = 0;

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

  std::string str1 = std::string("(3, 2)");
  std::string str2 = std::string("(2, 1)");
  std::string str3 = std::string("(4, 3)");
  std::string str4 = std::string("1");
  std::string str5 = std::string("9");

  CGarside::ComplexDualBraidParameter en(3, Index);
  CGarside::ComplexDualBraidUnderlying cf1(CGarside::ComplexDualBraidParameter(3, Index));
  CGarside::ComplexDualBraidUnderlying cf2(CGarside::ComplexDualBraidParameter(3, Index));
  cf1.OfString(str1);
  cf2.OfString(str2);
  CGarside::ComplexDualBraidUnderlying cf3 = cf1.Product(cf2);
  std::cout << "cf1: ";
  cf1.Debug(std::cout);
  std::cout << std::endl;
  std::cout << "cf2: ";
  cf2.Debug(std::cout);
  std::cout << std::endl;
  std::cout << "cf3: ";
  cf3.Debug(std::cout);
  std::cout << std::endl;
  int cfpart[3 * Index];
  cf3.AssignPartition(cfpart);
  std::cout << "cfpart: ";
  for (int i = 0; i <= 3 * Index; i++) {
    std::cout << cfpart[i] << " ";
  }
  std::cout << std::endl;
  cf3.OfPartition(cfpart);
  std::cout << "cf3: ";
  cf3.Debug(std::cout);
  std::cout << std::endl;
  std::string("");

  //CBraid::ArtinBraid cb = ToCBraid(b);
//
  double itime, ftime, exec_time;
//
  //itime = omp_get_wtime();
//
  //std::list<std::list<CBraid::ArtinBraid>> cscs = Braiding::SC(cb);
  
  //ftime = omp_get_wtime();

  //exec_time = ftime - itime;
//
  //std::cout << "Computing (Braiding) the SCS took " << exec_time * 1000 << " ms." << std::endl;
  //std::cout << "It has " << cscs.size() << " orbits, and " << cscs.front().size() * cscs.size()  << " elements." << std::endl;

  //t.Start();

  //SSS::SuperSummitSet<CGarside::ArtinBraid> sss = SSS::SSS(b);

  //t.Stop();

  //std::cout << "Computing the SSS, of size " << sss.Card() << ", took " << t.IntervalSec() * 1000 << " ms." << std::endl;

  itime = omp_get_wtime();

  SC::SlidingCircuitSet<CGarside::ArtinBraid> scs = SC::SCS(b);

  ftime = omp_get_wtime();

  exec_time = ftime - itime;

  //scs.Print(std::cout);

  std::cout << "Computing the SCS took " << exec_time * 1000 << " ms." << std::endl;
  std::cout << "It has " << scs.NumberOfOrbits() << " orbits, and " << scs.Card()  << " elements." << std::endl;

  return 0; 
}
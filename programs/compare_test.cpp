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

int Index = 8;
int CLength = 16;
int Count = 1;
int RandomSeed = 0;

CBraid::ArtinFactor ToCBraidFactor(const CGarside::ArtinBraidFactor &f)
{
  CBraid::ArtinFactor cf(f.GetParameter());
  CGarside::ArtinBraidUnderlying u = f.GetUnderlying();
  for (int i = 1; i <= f.GetParameter(); i++)
  {
    cf.At(i) = u.At(i);
  }
  return cf;
}

CBraid::ArtinBraid ToCBraid(const CGarside::ArtinBraid &b)
{
  CBraid::ArtinBraid cb(b.GetParameter());
  cb.LeftDelta = b.Inf();
  for (CGarside::ArtinBraid::ConstFactorItr it = b.FactorList.begin(); it != b.FactorList.end(); it++)
  {
    cb.FactorList.push_back(ToCBraidFactor(*it));
  }
  return cb;
}

int main()
{
  srand(RandomSeed);

  CGarside::ArtinBraid b(Index);

  b.Randomize(CLength);
  // b.Normalize();

  std::string str1 = std::string("(5, 2)");
  std::string str2 = std::string("(6, 5)");
  std::string str3 = std::string("(4, 3)");
  std::string str4 = std::string("1");
  std::string str5 = std::string("5");

  // int *cfpart = new int[3 * Index + 1];
  CGarside::ComplexDualBraidParameter en(3, Index);
  // CGarside::BandBraidFactor foo(Index);
  CGarside::ComplexDualBraid cfoo(en);
  // cfoo.Randomize(CLength);
  // CGarside::ComplexDualBraidUnderlying cf1(CGarside::ComplexDualBraidParameter(3, Index));
  // CGarside::ComplexDualBraidUnderlying cf2(CGarside::ComplexDualBraidParameter(3, Index));
  // cf1.OfString(str4);
  // cf1.AssignPartition(cfpart);
  // cf2.OfString(str5);
  // CGarside::ComplexDualBraidUnderlying cf3 = cf1.Product(cf2);
  // cf1.Delta();
  // cf3.OfString(str5);
  // std::cout << "cf1: ";
  // cf1.Debug(std::cout);
  // std::cout << std::endl;
  // std::cout << "cf1: " << cf1 << std::endl;
  // std::cout << "cf2: ";
  // cf2.Debug(std::cout);
  // std::cout << std::endl;
  // std::cout << "cf2: " << cf2 << std::endl;
  // std::cout << "cf3: ";
  // cf3.RightComplement(cf1).Debug(std::cout);
  // std::cout << std::endl;
  // std::cout << "cf3.RighComplement: ";
  // cf1.Delta();
  // cf2.OfString(str3);
  // cf3.RightComplement(cf1).LeftMeet(cf2).Debug(std::cout);
  // std::cout << std::endl;
  // cf3.RightComplement(cf1).AssignPartition(cfpart);
  // cf3.OfString(str5);
  // cf3.Debug(std::cout);
  // std::cout << std::endl;
  // cf3.AssignPartition(cfpart);
  // std::cout << "cfpart: ";

  // for (int i = 0; i <= 3 * Index; i++)
  // {
  // std::cout << cfpart[i] << " ";
  // }
  // std::cout << std::endl;
  // cf3.OfPartition(cfpart);
  // std::cout << "cf3 (of part): ";
  // cf3.Debug(std::cout);
  // std::cout << std::endl;
  // std::cout << cf3 << std::endl;
  //
  // cf1.OfString(str1);
  // cf2.OfString(str2);
  // cf3 = cf1.Product(cf2);
  // cf2.OfString(str5);
  // cf3 = cf3.Product(cf2);
  // cf3.Print(std::cout);
  // std::cout << std::endl;
  // cf3.Debug(std::cout);
  // std::cout << std::endl;
  // cf3.AssignPartition(cfpart);
  // std::cout << "cfpart: ";

  // for (int i = 0; i <= 3 * Index; i++)
  // {
  // std::cout << cfpart[i] << " ";
  // }
  // cf1.Delta();
  // cf2.OfString(str3);
  // cf2.Debug(std::cout);
  // std::cout << std::endl;
  // cf3.RightComplement(cf1).Debug(std::cout);
  // std::cout << std::endl;
  // cf3.DeltaConjugate(13);
  // cf3.Debug(std::cout);
  //
  // std::cout << std::endl;
  CGarside::ComplexDualBraidFactor cof(en);
  //

  std::vector<CGarside::ComplexDualBraidFactor> atoms = cof.Atoms();
  CGarside::ComplexDualBraid cob(en);

  for (int j = 0; j < Count; j++)
  {  for (int i = 0; i < CLength; i++)
  {
    cof = atoms[rand() % atoms.size()];
    cob.RightProduct(cof);
  }

  std::cout << cob << std::endl;
  cob.Debug(std::cout);
  std::cout << std::endl;

  SC::SCS(cob).Print(std::cout);

  }
  // delete[] cfpart;

  // CBraid::ArtinBraid cb = ToCBraid(b);
  //
  double itime, ftime, exec_time;
  //
  // itime = omp_get_wtime();
  //
  // std::list<std::list<CBraid::ArtinBraid>> cscs = Braiding::SC(cb);

  // ftime = omp_get_wtime();

  // exec_time = ftime - itime;
  //
  // std::cout << "Computing (Braiding) the SCS took " << exec_time * 1000 << " ms." << std::endl;
  // std::cout << "It has " << cscs.size() << " orbits, and " << cscs.front().size() * cscs.size()  << " elements." << std::endl;

  // t.Start();

  // SSS::SuperSummitSet<CGarside::ArtinBraid> sss = SSS::SSS(b);

  // t.Stop();

  // std::cout << "Computing the SSS, of size " << sss.Card() << ", took " << t.IntervalSec() * 1000 << " ms." << std::endl;

  // std::cout << b << std::endl;
  // b.Debug(std::cout);
  // std::cout << std::endl;

  itime = omp_get_wtime();

  //SC::SlidingCircuitSet<CGarside::ArtinBraid> scs = SC::SCS(b);

  ftime = omp_get_wtime();

  exec_time = ftime - itime;

  // scs.Print(std::cout);

  // std::cout << "Computing the SCS took " << exec_time * 1000 << " ms." << std::endl;
  // std::cout << "It has " << scs.NumberOfOrbits() << " orbits, and " << scs.Card() << " elements." << std::endl;

  return 0;
}
#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

#include "artin_braid.h"
#include "band_braid.h"
#include "garsidesc.h"
#include "cbraid.h"
#include "braiding.h"

#include "optarg.h"
#include "timecounter.h"

int Index = 8;
int CLength = 2;
int Count = 1;
int RandomSeed;

void PrintArray(int* array, int len) {
  for(int i = 1; i <= len; i++) {
    std::cout << array[i];
  } 
}

void DebugBallotSequence(char* array, int len) {
  int s = 0;
  for(int i = 1; i <= len; i++) {
    s += int(array[i]);
    std::cout << int(array[i]) << ", " << s << "; ";
  }
  std::cout << std::endl;
}

void TestBallotSequence(int n) {
  char ballot[2 * n + 1];
  cln::cl_I i;
  for (i = 1; i <= CBraid::GetCatalanNumber(n); i++) {
    CBraid::BallotSequence(n, i, ballot);
    DebugBallotSequence(ballot, 2 * n);
  }
}

void TestRandomBallotSequence(int n) {
  char ballot[2 * n + 1];
  cln::cl_I i = cln::random_I(cln::default_random_state, CBraid::GetCatalanNumber(Index)) + 1;
  CBraid::BallotSequence(n, i, ballot);
  DebugBallotSequence(ballot, 2 * n);
}

int main()
{ 
  char ballot[Index + 1];
  int revarray[Index + 1];
  int array[Index + 1];

  CBraid::BandPresentation pres(Index);

  TimeCounter t;

  CBraid::BandBraid cb(Index);

  CGarside::BandBraid cg(Index);

  CGarside::ArtinBraidFactor s1 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s2 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s3 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));

  CBraid::ArtinFactor cs1 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs2 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs3 = CBraid::ArtinFactor(Index, 0);

  CGarside::ArtinBraidFactor foo = s1, bar = s1;

  std::string str1 = std::string("1");
  s1.OfString(str1);
  std::string str2 = std::string("2");
  s2.OfString(str2);
  std::string str3 = std::string("3");
  s3.OfString(str3);

  str1 = std::string("1");
  str2 = std::string("2");
  str3 = std::string("3");
  cs1.OfString(str1);
  cs2.OfString(str2);
  cs3.OfString(str3);

  foo.Randomize();
  bar.Randomize();

  CGarside::BandBraidFactor a13 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a12 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a23 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));

  CGarside::BandBraidUnderlying ua13 = CGarside::BandBraidUnderlying(Index);

  CGarside::ArtinBraid cag(Index);

  CGarside::ArtinBraid scb1(s1);
  CGarside::ArtinBraid scb2(s2);
  CGarside::ArtinBraid scb3(s3);

  CBraid::ArtinBraid cscb1(cs1);
  CBraid::ArtinBraid cscb2(cs2);
  CBraid::ArtinBraid cscb3(cs3);
  CGarside::ArtinBraid scb = scb1 * scb2 * scb3 * scb1 * scb1 * scb2;
  scb.Normalize();
  CBraid::ArtinBraid cscb = cscb1 * cscb2 * cscb3 * cscb1 * cscb1 * cscb2;
  cscb.MakeLCF();

  SC::SlidingCircuitSet<CGarside::ArtinBraid> sc = SC::SCS(scb);
  sc.Print(std::cout);

  //std::string str13 = std::string("(3, 1)");
  //std::string str12 = std::string("(2, 1)");
  //std::string str23 = std::string("(3, 2)");

  s1.Identity();

  //a13.OfString(str13);
  //a12.OfString(str12);
  //a23.OfString(str23);

  a13.Identity();

  cb.Randomize(CLength);

  t.Start();

  Braiding::SC(cscb);

  t.Stop();

  std::cout << "Cbraid: Execution time: " << 1000 * t.IntervalSec() / Count << " ms." << std::endl;

  t.Start();

  //for (int i = 0; i < Count; i++)
  //{
  //  cag.Randomize(s1, CLength);
  //  cag.Normalize();
  //}

  t.Stop();

  cag.Randomize(CLength);

  t.Start();

  SC::SCS(scb);

  t.Stop();

  std::cout << "CGarside (Artin structure): Execution time: " << 1000 * t.IntervalSec() / Count << " ms." << std::endl;

  CGarside::MadeLeftWeighted = 0;

  t.Start();

  //for (int i = 0; i < Count; i++)
  //{
  //  cg.FactorList.clear(); // Potential memory leak? To be checked.
  //  cg.Delta = 0;
  //  cg.Randomize(a12, CLength);
  //  cg.Normalize();
  //}

  t.Stop();

    cg.Randomize(CLength);

    cg.Debug(std::cout);

    t.Start();
    
    cg.Normalize();

    t.Stop();

  std::cout << "CGarside (dual structure): Execution time: " << 1000 *  t.IntervalSec() / Count << " ms." << std::endl;

  std::cout << "CGarside (dual structure): " << CGarside::MadeLeftWeighted / Count << " calls to MakeLeftWeighted in average." << std::endl;

  return 0;
}
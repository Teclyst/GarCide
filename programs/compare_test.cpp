#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

#include "cbraid.h"
#include "artin_braid.h"
#include "band_braid.h"

#include "optarg.h"
#include "timecounter.h"

int Index = 50;
int CLength = 200;
int Count = 1;
int RandomSeed;

void PrintArray(int* array, int len) {
  for(int i = 1; i <= len; i++) {
    std::cout << array[i];
  } 
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

  CGarside::ArtinBraidFactor foo = s1, bar = s1;

  std::string str1 = std::string("1");

  foo.OfString(str1);

  foo.Print(std::cout);
  bar.Print(std::cout);
  std::cout << std::endl;

  while (true) {};

  CGarside::BandBraidFactor a13 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a12 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a23 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));

  CGarside::BandBraidUnderlying ua13 = CGarside::BandBraidUnderlying(Index);

  CGarside::ArtinBraid cag(Index);

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

  for (int i = 0; i < Count; i++)
  {
    cb.MakeLCF();
  }

  t.Stop();

  std::cout << "Cbraid: Execution time: " << 1000 * t.IntervalSec() / Count << " ms." << std::endl;

  std::cout << "Cbraid: " << CBraid::MadeLeftWeighted / Count << " calls to MakeLeftWeighted in average." << std::endl;

  t.Start();

  //for (int i = 0; i < Count; i++)
  //{
  //  cag.Randomize(s1, CLength);
  //  cag.Normalize();
  //}

  t.Stop();

    cag.Randomize(CLength);

    cag.Print(std::cout);

    t.Start();

    cag.Normalize();

  t.Stop();

  cag.Print(std::cout);

  std::cout << "CGarside (Artin structure): Execution time: " << 1000 *  t.IntervalSec() / Count << " ms." << std::endl;

  std::cout << "CGarside (Artin structure): " << CGarside::MadeLeftWeighted / Count << " calls to MakeLeftWeighted in average." << std::endl;

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

    cg.Print(std::cout);

    t.Start();
    
    cg.Normalize();

    t.Stop();

    cg.Print(std::cout);

  std::cout << "CGarside (dual structure): Execution time: " << 1000 *  t.IntervalSec() / Count << " ms." << std::endl;

  std::cout << "CGarside (dual structure): " << CGarside::MadeLeftWeighted / Count << " calls to MakeLeftWeighted in average." << std::endl;

  return 0;
}
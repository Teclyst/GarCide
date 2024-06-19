#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

#include "cbraid.h"
#include "artin_braid.h"
#include "band_braid.h"

#include "optarg.h"
#include "timecounter.h"

int Index = 150;
int CLength = 1000;
int Count = 1;
int RandomSeed;

int main()
{

  TimeCounter t;

  CBraid::BandBraid cb(Index);

  CGarside::BandBraid cg;

  CGarside::ArtinBraidFactor s1 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s2 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s3 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));

  CGarside::BandBraidFactor a13 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a12 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a23 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));

  CGarside::BandBraidUnderlying ua13 = CGarside::BandBraidUnderlying(Index);

  CGarside::ArtinBraid cag;

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

    cag.Randomize(s1, CLength);
    t.Start();

    cag.Normalize();

  t.Stop();

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

    cg.FactorList.clear(); // Potential memory leak? To be checked.
    cg.Delta = 0;
    cg.Randomize(a12, CLength);

    t.Start();
    
    cg.Normalize();

    t.Stop();

  std::cout << "CGarside (dual structure): Execution time: " << 1000 *  t.IntervalSec() / Count << " ms." << std::endl;

  std::cout << "CGarside (dual structure): " << CGarside::MadeLeftWeighted / Count << " calls to MakeLeftWeighted in average." << std::endl;

  return 0;
}
#include <string>
#include <iostream>
#include <ctime>
#include <cmath>

#include "cbraid.h"
#include "artin_braid.h"

#include "optarg.h"
#include "timecounter.h"

int Index = 3;
int CLength = 5;
int Count = 100;
int RandomSeed;

int main()
{

  TimeCounter t;

  CBraid::ArtinBraid cb(Index);

  CGarside::ArtinBraid cg;

  CGarside::ArtinBraidFactor foo = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));

  foo.Identity();

  t.Start();

  for (int i = 0; i < Count; i++)
  {
    cb.Randomize(CLength);
    cb.MakeLCF();
  }

  t.Stop();

  std::cout << "Cbraid: Execution time: " << t.IntervalSec() << " s." << std::endl;

  std::cout << "Cbraid: " << CBraid::MadeLeftWeighted << " calls to MakeLeftWeighted." << std::endl;

  t.Start();

  for (int i = 0; i < Count; i++)
  {
    cg.FactorList.clear(); // Potential memory leak? To be checked.
    cg.Delta = 0;
    for (int i = 0; i < CLength; i++)
    {
      CGarside::ArtinBraidFactor f = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
      f.Delta();
      f.Randomize();
      CGarside::ArtinBraidFactor(f).Print(std::cout);
      std::cout << " . ";
      cg.FactorList.push_back(f);
    }
    std::cout << std::endl;
    cg.Normalize();
    cg.Print(std::cout);
  }

  t.Stop();

  std::cout << "CGarside: Execution time: " << t.IntervalSec() << " s." << std::endl;

  std::cout << "CGarside: " << CGarside::MadeLeftWeighted << " calls to MakeLeftWeighted." << std::endl;

  return 0;
}
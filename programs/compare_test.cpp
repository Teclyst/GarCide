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

int Index = 5;
int CLength = 1;
int Count = 1;
int RandomSeed;

int main()
{ 
  char ballot[Index + 1];
  int revarray[Index + 1];
  int array[Index + 1];

  CBraid::BandPresentation pres(Index);

  TimeCounter t;

  CBraid::BandBraid cb(Index);

  CGarside::BandBraid cg(Index);

  /*
  CGarside::ArtinBraidFactor s1 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s2 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s3 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s4 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s5 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s6 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));
  CGarside::ArtinBraidFactor s7 = CGarside::ArtinBraidFactor(CGarside::ArtinBraidUnderlying(Index));


  CBraid::ArtinFactor cs1 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs2 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs3 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs4 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs5 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs6 = CBraid::ArtinFactor(Index, 0);
  CBraid::ArtinFactor cs7 = CBraid::ArtinFactor(Index, 0);

  CGarside::ArtinBraidFactor foo = s1, bar = s1;

  std::string str1 = std::string("1");
  s1.OfString(str1);
  std::string str2 = std::string("2");
  s2.OfString(str2);
  std::string str3 = std::string("3");
  s3.OfString(str3);
  std::string str4 = std::string("4");
  s4.OfString(str4);
  std::string str5 = std::string("5");
  s5.OfString(str5);
  std::string str6 = std::string("6");
  s6.OfString(str6);
  std::string str7 = std::string("7");
  s7.OfString(str7);

  str1 = std::string("1");
  cs1.OfString(str1);
  str2 = std::string("2");
  cs2.OfString(str2);
  str3 = std::string("3");
  cs3.OfString(str3);
  str4 = std::string("4");
  cs4.OfString(str4);
  str5 = std::string("5");
  cs5.OfString(str5);
  str6 = std::string("6");
  cs6.OfString(str6);
  str7 = std::string("7");
  cs7.OfString(str7);

  CGarside::BandBraidFactor a13 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a12 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));
  CGarside::BandBraidFactor a23 = CGarside::BandBraidFactor(CGarside::BandBraidUnderlying(Index));

  CGarside::BandBraidUnderlying ua13 = CGarside::BandBraidUnderlying(Index);

  CGarside::ArtinBraid cag(Index);

  CGarside::ArtinBraid scb1(s1);
  CGarside::ArtinBraid scb2(s2);
  CGarside::ArtinBraid scb3(s3);
  CGarside::ArtinBraid scb4(s4);
  CGarside::ArtinBraid scb5(s5);
  CGarside::ArtinBraid scb6(s6);
  CGarside::ArtinBraid scb7(s7);

  CBraid::ArtinBraid cscb1(cs1);
  CBraid::ArtinBraid cscb2(cs2);
  CBraid::ArtinBraid cscb3(cs3);
  CBraid::ArtinBraid cscb4(cs4);
  CBraid::ArtinBraid cscb5(cs5);
  CBraid::ArtinBraid cscb6(cs6);
  CBraid::ArtinBraid cscb7(cs7);
  CGarside::ArtinBraid scb = scb7 * scb6 * scb5 * scb4 * scb3 * scb2 * scb1;
  scb.Normalize();
  CBraid::ArtinBraid cscb = cscb7 * cscb6 * cscb5 * cscb4 * cscb3 * cscb2 * cscb1;
  cscb.MakeLCF();

  //std::string str13 = std::string("(3, 1)");
  //std::string str12 = std::string("(2, 1)");
  //std::string str23 = std::string("(3, 2)");

  t.Start();

  for (int i = 0; i < Count; i++) {

  Braiding::SC(cscb);

  }
  t.Stop();

  std::cout << "Cbraid: Execution time: " << 1000 * t.IntervalSec() / Count << " ms." << std::endl;

  t.Start();

  //for (int i = 0; i < Count; i++)
  //{
  //  cag.Randomize(s1, CLength);
  //  cag.Normalize();
  //}

  t.Stop();

  t.Start();

  for (int i = 0; i < Count; i++) {

  SC::SCS(scb);

  }

  t.Stop();

  std::cout << "CGarside (Artin structure): Execution time: " << 1000 * t.IntervalSec() / Count << " ms." << std::endl;

  */
  //cg.Randomize(CLength);
  //cg.Normalize();
  //cg.Print(std::cout);
  //SC::SCS(cg).Print(std::cout);

  CGarside::BDualBraidFactor s1 = CGarside::BDualBraidFactor(Index);
  CGarside::BDualBraidFactor s2 = CGarside::BDualBraidFactor(Index);
  std::string str1 = std::string("1");
  s1.OfString(str1);
  std::string str2 = std::string("(2, 1)");
  s2.OfString(str2);
  (s2 * s1).Debug(std::cout);
  (s2 * s1).Print(std::cout);
  SSS::SSS(CGarside::Braid(s1)).Debug(std::cout);
  USS::USS(CGarside::Braid(s1)).Debug(std::cout);

  CGarside::IDualBraidFactor is1 = CGarside::IDualBraidFactor(Index);
  CGarside::IDualBraidFactor is2 = CGarside::IDualBraidFactor(Index);
  std::string istr1 = std::string("0");
  is1.OfString(istr1);
  std::string istr2 = std::string("1");
  is2.OfString(istr2);
  (is2).Debug(std::cout);
  (is1).Debug(std::cout);
  std::cout << std::endl;
  SSS::SSS(CGarside::Braid(is1)).Debug(std::cout);

  return 0;
}
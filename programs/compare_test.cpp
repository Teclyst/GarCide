#include <cmath>
#include <ctime>
#include <iostream>
#include <string>

#include "artin_braid.h"
#include "band_braid.h"
#include "garsidesc.h"

int Index = 8;
int CLength = 8;
int Count = 1;
int RandomSeed = 0;

Utility::IndentedOStream ind_cout(std::cout);

// CBraid::ArtinFactor ToCBraidFactor(const CGarside::ArtinBraidFactor &f)
//{
//   CBraid::ArtinFactor cf(f.GetParameter());
//   CGarside::ArtinBraidUnderlying u = f.GetUnderlying();
//   for (int i = 1; i <= f.GetParameter(); i++)
//   {
//     cf.At(i) = u.At(i);
//   }
//   return cf;
// }
//
// CBraid::ArtinBraid ToCBraid(const CGarside::ArtinBraid &b)
//{
//   CBraid::ArtinBraid cb(b.GetParameter());
//   cb.LeftDelta = b.Inf();
//   for (CGarside::ArtinBraid::ConstFactorItr it = b.FactorList.begin(); it !=
//   b.FactorList.end(); it++)
//   {
//     cb.FactorList.push_back(ToCBraidFactor(*it));
//   }
//   return cb;
// }

int main() {
    std::string str = std::string("D ^ 5 . 2 4 3 2 1 2  3 2 1 2 1  2 2 1  2 2 "
                                  "1 2 1 7 5 7 5 4 3 4 5 6 6 5 4 6");

    CGarside::ArtinBraid b(Index);

    try {
        b.OfString(str);
    } catch (Utility::InvalidStringError inval) {
        std::cout << inval.error_source << std::endl;
    }

    SC::SCS(b).Print(ind_cout);
    ind_cout << Utility::EndLine();
    return 0;
}
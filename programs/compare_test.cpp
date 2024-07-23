// #include <cmath>
// #include <ctime>
#include <iostream>
#include <string>

#include "artin_braid.h"
// #include "band_braid.h"
#include "braiding.h"
#include "cbraid.h"
// #include "dual_complex_reflection.h"
// #include "standard_complex_reflection.h"
#include "garsidesc.h"
#include "garsideuss.h"
// #include "octaedral_braid.h"

int Index = 8;
int E = 3;
int CLength = 4;
int Count = 20;
int RandomSeed = 0;

// CBraid::ArtinFactor ToCBraidFactor(const cgarside::ArtinBraidFactor &f) {
    // CBraid::ArtinFactor cf(f.GetParameter());
    // cgarside::ArtinBraidUnderlying u = f.GetUnderlying();
    // for (int i = 1; i <= f.GetParameter(); i++) {
        // cf.At(i) = u.At(i);
    // }
    // return cf;
// }
// 
// CBraid::ArtinBraid ToCBraid(const cgarside::ArtinBraid &b) {
    // CBraid::ArtinBraid cb(b.GetParameter());
    // cb.LeftDelta = b.Inf();
    // for (cgarside::ArtinBraid::ConstFactorItr it = b.FactorList.begin();
        //  it != b.FactorList.end(); it++) {
        // cb.FactorList.push_back(ToCBraidFactor(*it));
    // }
    // return cb;
// }
// 
// cgarside::ArtinBraidFactor OfCBraidFactor(const CBraid::ArtinFactor &f) {
    // cgarside::ArtinBraidUnderlying u(f.Index());
    // for (int i = 1; i <= f.Index(); i++) {
        // u.At(i) = f.At(i);
    // }
    // return cgarside::ArtinBraidFactor(u);
// }
// 
// cgarside::ArtinBraid OfCBraid(const CBraid::ArtinBraid &b) {
    // cgarside::ArtinBraid cb(b.Index());
    // cb.Delta = b.LeftDelta;
    // for (CBraid::ArtinBraid::ConstFactorItr it = b.FactorList.begin();
        //  it != b.FactorList.end(); it++) {
        // cb.RightProduct(OfCBraidFactor(*it));
    // }
    // return cb;
// }
// 
// void PrintFactorList(const std::list<cgarside::ArtinBraidFactor> &l) {
    // for (typename std::list<cgarside::ArtinBraidFactor>::const_iterator it =
            //  l.begin();
        //  it != l.end(); it++) {
        // (*it).Print(cgarside::ind_cout);
        // cgarside::ind_cout << cgarside::EndLine();
    // }
// }
// 
// void PrintCBraidFactorList(const std::list<CBraid::ArtinFactor> &l) {
    // for (typename std::list<CBraid::ArtinFactor>::const_iterator it = l.begin();
        //  it != l.end(); it++) {
        // OfCBraidFactor(*it).Print(cgarside::ind_cout);
        // cgarside::ind_cout << cgarside::EndLine();
    // }
// }
// 
// void PrintList(const std::list<cgarside::ArtinBraid> &l) {
    // for (typename std::list<cgarside::ArtinBraid>::const_iterator it =
            //  l.begin();
        //  it != l.end(); it++) {
        // (*it).Print(cgarside::ind_cout);
        // cgarside::ind_cout << cgarside::EndLine();
    // }
// }
// 
// void CBraidPrintList(const std::list<CBraid::ArtinBraid> &l) {
    // for (typename std::list<CBraid::ArtinBraid>::const_iterator it = l.begin();
        //  it != l.end(); it++) {
        // OfCBraid(*it).Print(cgarside::ind_cout);
        // cgarside::ind_cout << cgarside::EndLine();
    // }
// }

int main() {
    std::string str = std::string(
        "1 2 1 3 2 1 4 5 4 6 5 4 3 2 1 7 6 5 4 3 2 1 . 1 2 1 3 2 1 4 3 2 1 5 4 "
        "3 2 6 5 4 3 2 7 6 5 4 3 2 1 . 3 4 3 5 4 3 . 4 5 4 3 6 5");

    // cgarside::ComplexDualBraid b(cgarside::ComplexDualBraidParameter(E,
    // Index)); cgarside::ComplexStandardBraid
    // c(cgarside::ComplexStandardBraidParameter(Index, E));
    // cgarside::ComplexStandardBraidFactor
    // f(cgarside::ComplexStandardBraidParameter(Index, E));
    // cgarside::ComplexStandardBraidFactor
    // g(cgarside::ComplexStandardBraidParameter(Index, E));

    cgarside::artin::ArtinBraid b(Index);
    cgarside::artin::ArtinBraid b_rcf(Index);

    cgarside::artin::ArtinFactor f(Index);

    std::vector<cgarside::artin::ArtinBraid> t, t_rcf;

    try {
        b.OfString(str);
    } catch (cgarside::InvalidStringError inval) {
        std::cout << inval.error_source << std::endl;
    }

    cgarside::sliding_circuit::SCS(b).print(cgarside::ind_cout);
    cgarside::ind_cout << "b is "<< cgarside::artin::thurston_type(b) << "." << cgarside::EndLine();
    return 0;
}
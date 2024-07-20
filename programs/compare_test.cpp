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
int Count = 1000;
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
    for (CGarside::ArtinBraid::ConstFactorItr it = b.FactorList.begin();
         it != b.FactorList.end(); it++) {
        cb.FactorList.push_back(ToCBraidFactor(*it));
    }
    return cb;
}

CGarside::ArtinBraidFactor OfCBraidFactor(const CBraid::ArtinFactor &f) {
    CGarside::ArtinBraidUnderlying u(f.Index());
    for (int i = 1; i <= f.Index(); i++) {
        u.At(i) = f.At(i);
    }
    return CGarside::ArtinBraidFactor(u);
}

CGarside::ArtinBraid OfCBraid(const CBraid::ArtinBraid &b) {
    CGarside::ArtinBraid cb(b.Index());
    cb.Delta = b.LeftDelta;
    for (CBraid::ArtinBraid::ConstFactorItr it = b.FactorList.begin();
         it != b.FactorList.end(); it++) {
        cb.RightProduct(OfCBraidFactor(*it));
    }
    return cb;
}

int main() {
    std::string str = std::string(
        "1 2 1 3 2 1 4 5 4 6 5 4 3 2 1 7 6 5 4 3 2 1 . 1 2 1 3 2 1 4 3 2 1 5 4 "
        "3 2 6 5 4 3 2 7 6 5 4 3 2 1 . 3 4 3 5 4 3 . 4 5 4 3 6 5");

    // CGarside::ComplexDualBraid b(CGarside::ComplexDualBraidParameter(E,
    // Index)); CGarside::ComplexStandardBraid
    // c(CGarside::ComplexStandardBraidParameter(Index, E));
    // CGarside::ComplexStandardBraidFactor
    // f(CGarside::ComplexStandardBraidParameter(Index, E));
    // CGarside::ComplexStandardBraidFactor
    // g(CGarside::ComplexStandardBraidParameter(Index, E));

    CGarside::ArtinBraid b(Index);
    CGarside::ArtinBraid b_rcf(Index);

    CGarside::ArtinBraidFactor f(Index);

    std::vector<CGarside::ArtinBraid> t, t_rcf;

    for (int i = 0; i < Count; i++) {
        b.Randomize(CLength);
        f.Randomize();
        b.Normalize();
        b = USS::SendToUSS(b);
        b_rcf = b;
        b_rcf.MakeRCFFromLCF();
        if (!f.IsIdentity() && USS::MainPullback(b, b_rcf, f) !=
                                   OfCBraidFactor(Braiding::MainPullback(
                                       ToCBraid(b), ToCBraidFactor(f)))) {
            b.Print(CGarside::ind_cout);
            CGarside::ind_cout << CGarside::EndLine();
            USS::MainPullback(b, b_rcf, f).Print(CGarside::ind_cout);
            CGarside::ind_cout << CGarside::EndLine();
            OfCBraidFactor(
                Braiding::MainPullback(ToCBraid(b), ToCBraidFactor(f)))
                .Print(CGarside::ind_cout);
            CGarside::ind_cout << CGarside::EndLine();
        }
    }

    CGarside::ind_cout << "Passed all tests!" << CGarside::EndLine();

    try {
        b.OfString(str);
    } catch (CGarside::InvalidStringError inval) {
        std::cout << inval.error_source << std::endl;
    }

    SSS::SendToSSS(b);
    CGarside::ind_cout << CGarside::EndLine();

    Braiding::SendToSSS(ToCBraid(b));
    CGarside::ind_cout << CGarside::EndLine();
    return 0;
}
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>

#include "artin_braid.h"
#include "band_braid.h"
#include "braiding.h"
#include "cbraid.h"
#include "dual_complex_reflection.h"
#include "standard_complex_reflection.h"
#include "garsidesc.h"
#include "octaedral_braid.h"

int Index = 8;
int E = 3;
int CLength = 8;
int Count = 1000;
int RandomSeed = 0;

CBraid::BandFactor ToCBraidFactor(const CGarside::BandBraidFactor &f) {
    CBraid::BandFactor cf(f.GetParameter());
    CGarside::BandBraidUnderlying u = f.GetUnderlying();
    for (int i = 1; i <= f.GetParameter(); i++) {
        cf.At(i) = u.At(i);
    }
    return cf;
}

CBraid::BandBraid ToCBraid(const CGarside::BandBraid &b) {
    CBraid::BandBraid cb(b.GetParameter());
    cb.LeftDelta = b.Inf();
    for (CGarside::BandBraid::ConstFactorItr it = b.FactorList.begin();
         it != b.FactorList.end(); it++) {
        cb.FactorList.push_back(ToCBraidFactor(*it));
    }
    return cb;
}

CGarside::BandBraidFactor OfCBraidFactor(const CBraid::BandFactor &f) {
    CGarside::BandBraidUnderlying u(f.Index());
    for (int i = 1; i <= f.Index(); i++) {
        u.At(i) = f.At(i);
    }
    return CGarside::BandBraidFactor(u);
}

CGarside::BandBraid OfCBraid(const CBraid::BandBraid &b) {
    CGarside::BandBraid cb(b.Index());
    cb.Delta = b.LeftDelta;
    for (CBraid::BandBraid::ConstFactorItr it = b.FactorList.begin();
         it != b.FactorList.end(); it++) {
        cb.RightProduct(OfCBraidFactor(*it));
    }
    return cb;
}

int main() {
    std::string str =
        std::string("D ^ 5 s4 s _ 3 s6 t6 ");

    CGarside::ComplexStandardBraid b(CGarside::ComplexStandardBraidParameter(E, Index));
    CGarside::ComplexStandardBraid c(CGarside::ComplexStandardBraidParameter(Index, E));
    CGarside::ComplexStandardBraidFactor f(CGarside::ComplexStandardBraidParameter(Index, E));
    CGarside::ComplexStandardBraidFactor g(CGarside::ComplexStandardBraidParameter(Index, E));

    for (int i = 0; i < Count; i++) {
        //b.Randomize(CLength);
        //b.Normalize();
        // c.Randomize(CLength);
        // c.Normalize();
        // if (b.LeftJoin(c) !=
        //     OfCBraid(Braiding::LeftWedge(ToCBraid(b), ToCBraid(c)))) {
        //     CGarside::ind_cout << "Error!";
        // }
    }

    try {
        b.OfString(str);
    } catch (CGarside::InvalidStringError inval) {
        std::cout << inval.error_source << std::endl;
    }

    SC::SCS(b).Debug(CGarside::ind_cout);
    CGarside::ind_cout << CGarside::EndLine();
    CGarside::ind_cout << CGarside::EndLine();
    CGarside::ind_cout << b << CGarside::EndLine();
    b.Debug(CGarside::ind_cout);
    CGarside::ind_cout << CGarside::EndLine();
    return 0;
}
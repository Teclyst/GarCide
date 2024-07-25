#include "standard_complex_reflection.h"
#include <iostream>
#include <string>

namespace cgarside {

namespace standard_complex {

void Parameter::Print(IndentedOStream &os) const {
    os << "(e: " << e << ", n: " << n << ")";
}

void Underlying::Debug(IndentedOStream &os) const {
    os << "{   ";
    os.Indent(4);
    os << "PresentationParameter:";
    os.Indent(4);
    os << EndLine() << GetParameter();
    os.Indent(-4);
    os << EndLine();
    os << "PermutationTable:";
    os.Indent(4);
    os << EndLine();
    os << "[";
    for (sint16 i = 0; i < GetParameter().n; i++) {
        os << PermutationTable[i] << ", ";
    }
    os << PermutationTable[GetParameter().n];
    os << "]";
    os.Indent(-4);
    os << EndLine();
    os << "CoefficientTable:";
    os.Indent(4);
    os << EndLine();
    os << "[";
    for (sint16 i = 0; i < GetParameter().n; i++) {
        os << CoefficientTable[i] << ", ";
    }
    os << CoefficientTable[GetParameter().n];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

Parameter Underlying::GetParameter() const { return PresentationParameter; }

Underlying::ParameterType
Underlying::parameter_of_string(const std::string &str) {
    std::smatch match;

    if (std::regex_match(str, match,
                         std::regex{"[\\s\\t]*\\([\\s\\t]*(" + number_regex +
                                    ")[\\s\\t]*,?[\\s\\t]*(" + number_regex +
                                    ")[\\s\\t]*\\)[\\s\\t]*"},
                         std::regex_constants::match_continuous)) {
        sint16 e, n;
        try {
            e = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Parameter e is too big!\n" +
                                     match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        try {
            n = std::stoi(match[2]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Parameter n is too big!\n" +
                                     match.str(2) +
                                     " can not be converted to a C++ integer.");
        }
        if ((2 <= e) && (2 <= n) && (n <= MaxN)) {
            return Parameter(e, n);
        } else if (2 > e) {
            throw InvalidStringError("e should be at least 2!");
        } else if (2 > n) {
            throw InvalidStringError("n should be at least 2!");
        } else {
            throw InvalidStringError("n is too big!\n" + match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MaxN) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

sint16 Underlying::LatticeHeight() const {
    sint16 n = GetParameter().n;
    return n * (n + 1);
}

Underlying::Underlying(Parameter p)
    : PresentationParameter(p), PermutationTable(p.n), CoefficientTable(p.n) {}

void Underlying::Print(IndentedOStream &os) const {
    thread_local sint16 dir_perm[MaxN];
    Underlying copy = *this;
    copy.Direct(dir_perm);
    sint16 n = GetParameter().n, e = GetParameter().e;
    bool is_first = true;
    for (sint16 i = 2; i <= n; i++) {
        for (sint16 j = i; j > 2; j--) {
            if (copy.IsSLeftDivisor(dir_perm, j)) {
                if (is_first) {
                    is_first = false;
                } else {
                    os << " ";
                }
                copy.SLeftMultiply(dir_perm, j);
                os << "s" << j;
            } else {
                goto end_subword;
            }
        }

        if (copy.IsTLeftDivisor(dir_perm, 1)) {
            if (is_first) {
                is_first = false;
            } else {
                os << " ";
            }
            copy.TLeftMultiply(dir_perm, 1);
            os << "t" << 1;

            if (copy.IsTLeftDivisor(dir_perm, 0)) {
                copy.TLeftMultiply(dir_perm, 0);
                os << " t" << 0;

                for (sint16 j = 3; j < i + 1; j++) {
                    if (copy.IsSLeftDivisor(dir_perm, j)) {
                        copy.SLeftMultiply(dir_perm, j);
                        os << " s" << j;
                    } else {
                        goto end_subword;
                    }
                }
            }
            goto end_subword;
        }

        for (sint16 k = 0; k < e; k++) {
            if (copy.IsTLeftDivisor(dir_perm, k)) {
                if (is_first) {
                    is_first = false;
                } else {
                    os << " ";
                }
                copy.TLeftMultiply(dir_perm, k);
                os << "t" << k;
                goto end_subword;
            }
        }
    end_subword:;
    }
}

void Underlying::OfString(const std::string &str, size_t &pos) {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::smatch match;
    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        Delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"([st])[\\s\\t]*_?[\\s\\t]*(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[2]);
        pos += match[0].length();
        if (match[1] == "s") {
            if ((3 <= i) && (i <= n)) {
                // s_i.
                Identity();
                PermutationTable[i - 2] = i - 1;
                PermutationTable[i - 1] = i - 2;
            } else {
                throw InvalidStringError(
                    "Invalid index for s type generator!\n" +
                    std::to_string(i) + " is not in [3, " + std::to_string(n) +
                    "].");
            }
        } else {
            // t_i.
            i = Rem(i, e);
            Identity();
            PermutationTable[0] = 1;
            PermutationTable[1] = 0;
            CoefficientTable[1] = i;
            CoefficientTable[0] = e - ((i == 0) ? e : i);
        }
    } else {
        throw InvalidStringError(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('s' | 't') '_'? "
            "Z | 'D',\n where Z matches integers, and ignoring whitespaces.");
        ;
    }
}

void Underlying::Direct(sint16 *dir_perm) const {
    for (sint16 i = 0; i < GetParameter().n; i++) {
        dir_perm[PermutationTable[i]] = i;
    }
};

Underlying Underlying::LeftMeet(const Underlying &b) const {
    thread_local sint16 dir_perm_a[MaxN], dir_perm_b[MaxN], dir_perm_meet[MaxN];
    Underlying a_copy = *this;
    Underlying b_copy = b;
    Underlying meet(GetParameter());
    meet.Identity();
    a_copy.Direct(dir_perm_a);
    b_copy.Direct(dir_perm_b);
    meet.Direct(dir_perm_meet);
    sint16 n = GetParameter().n, e = GetParameter().e;
    for (sint16 i = 2; i <= n; i++) {
        for (sint16 j = i; j > 2; j--) {
            if (a_copy.IsSLeftDivisor(dir_perm_a, j) &&
                b_copy.IsSLeftDivisor(dir_perm_b, j)) {
                meet.SRightMultiply(dir_perm_meet, j);
                a_copy.SLeftMultiply(dir_perm_a, j);
                b_copy.SLeftMultiply(dir_perm_b, j);
            } else {
                goto end_subword;
            }
        }

        if (a_copy.IsTLeftDivisor(dir_perm_a, 0) &&
            b_copy.IsTLeftDivisor(dir_perm_b, 0)) {

            meet.TRightMultiply(dir_perm_meet, 0);
            a_copy.TLeftMultiply(dir_perm_a, 0);
            b_copy.TLeftMultiply(dir_perm_b, 0);

            if (a_copy.IsTLeftDivisor(dir_perm_a, e - 1) &&
                b_copy.IsTLeftDivisor(dir_perm_b, e - 1)) {
                meet.TRightMultiply(dir_perm_meet, e - 1);
                a_copy.TLeftMultiply(dir_perm_a, e - 1);
                b_copy.TLeftMultiply(dir_perm_b, e - 1);

                for (sint16 j = 3; j < i + 1; j++) {
                    if (a_copy.IsSLeftDivisor(dir_perm_a, j) &&
                        b_copy.IsSLeftDivisor(dir_perm_b, j)) {
                        meet.SRightMultiply(dir_perm_meet, j);
                        a_copy.SLeftMultiply(dir_perm_a, j);
                        b_copy.SLeftMultiply(dir_perm_b, j);
                    } else {
                        goto end_subword;
                    }
                }
            }
            goto end_subword;
        }

        for (sint16 k = 1; k < e; k++) {
            if (a_copy.IsTLeftDivisor(dir_perm_a, k) &&
                b_copy.IsTLeftDivisor(dir_perm_b, k)) {
                meet.TRightMultiply(dir_perm_meet, k);
                a_copy.TLeftMultiply(dir_perm_a, k);
                b_copy.TLeftMultiply(dir_perm_b, k);
                goto end_subword;
            }
        }
    end_subword:;
    }
    return meet;
}

void Underlying::Identity() {
    for (sint16 i = 0; i < GetParameter().n; i++) {
        PermutationTable[i] = i;
        CoefficientTable[i] = 0;
    }
}

void Underlying::Delta() {
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 1; i < n; i++) {
        PermutationTable[i] = i;
        CoefficientTable[i] = 1;
    }
    PermutationTable[0] = 0;
    CoefficientTable[0] = Rem(-n + 1, e);
}

bool Underlying::Compare(const Underlying &b) const {
    sint16 i;
    for (i = 0; i < GetParameter().n; i++) {
        if ((PermutationTable[i] != b.PermutationTable[i]) ||
            (CoefficientTable[i] != b.CoefficientTable[i])) {
            return false;
        }
    }
    return true;
};

Underlying Underlying::Inverse() const {
    Underlying f = Underlying(GetParameter());
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 0; i < n; i++) {
        f.PermutationTable[PermutationTable[i]] = i;
        f.CoefficientTable[PermutationTable[i]] =
            CoefficientTable[i] == 0 ? 0 : e - CoefficientTable[i];
    }
    return f;
};

Underlying Underlying::Product(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 0; i < n; i++) {
        f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
        f.CoefficientTable[i] = Rem(
            b.CoefficientTable[PermutationTable[i]] + CoefficientTable[i], e);
    }
    return f;
};

Underlying Underlying::LeftComplement(const Underlying &b) const {
    return b.Product(Inverse());
};

Underlying Underlying::RightComplement(const Underlying &b) const {
    return Inverse().Product(b);
};

void Underlying::DeltaConjugate(sint16 k) {
    // Delta is diagonal, and acts almost homothetically.
    // Therefore conjugating by some power of it does nothing on most
    // coefficients.
    sint16 n = GetParameter().n, e = GetParameter().e;

    // In this case `*this` commutes with Delta.
    if (PermutationTable[0] == 0) {
        return;
    }

    // Otherwise the two non trivial coefficient are 0 and the i such that
    // `PermutationTable[i] == 0`.
    CoefficientTable[0] = Rem(CoefficientTable[0] + k * n, e);
    for (sint16 i = 1; i < n; i++) {
        if (PermutationTable[i] == 0) {
            CoefficientTable[i] = Rem(CoefficientTable[i] - k * n, e);
            return;
        }
    }
}

void Underlying::Randomize() { throw NonRandomizable(); }

std::vector<Underlying> Underlying::Atoms() const {
    Parameter p = GetParameter();
    std::vector<Underlying> atoms;
    Underlying atom = Underlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 2; i <= n - 1; i++) {
        // s_(i+1).
        atom.Identity();
        atom.PermutationTable[i - 1] = i;
        atom.PermutationTable[i] = i - 1;
        atoms.push_back(atom);
    }
    for (sint16 k = 0; k < e; k++) {
        // t_k.
        atom.Identity();
        atom.PermutationTable[0] = 1;
        atom.PermutationTable[1] = 0;
        atom.CoefficientTable[1] = k;
        atom.CoefficientTable[0] = e - ((k == 0) ? e : k);
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace standard_complex

template <>
IndentedOStream &IndentedOStream::operator<< <standard_complex::Parameter>(
    const standard_complex::Parameter &p) {
    p.Print(*this);
    return *this;
};

} // namespace cgarside
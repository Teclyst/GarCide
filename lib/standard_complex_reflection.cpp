#include "standard_complex_reflection.h"
#include <iostream>
#include <string>

namespace cgarside {

namespace standard_complex {

void Parameter::print(IndentedOStream &os) const {
    os << "(e: " << e << ", n: " << n << ")";
}

void Underlying::debug(IndentedOStream &os) const {
    os << "{   ";
    os.Indent(4);
    os << "PresentationParameter:";
    os.Indent(4);
    os << EndLine() << get_parameter();
    os.Indent(-4);
    os << EndLine();
    os << "permutation_table:";
    os.Indent(4);
    os << EndLine();
    os << "[";
    for (sint16 i = 0; i < get_parameter().n; i++) {
        os << permutation_table[i] << ", ";
    }
    os << permutation_table[get_parameter().n];
    os << "]";
    os.Indent(-4);
    os << EndLine();
    os << "CoefficientTable:";
    os.Indent(4);
    os << EndLine();
    os << "[";
    for (sint16 i = 0; i < get_parameter().n; i++) {
        os << CoefficientTable[i] << ", ";
    }
    os << CoefficientTable[get_parameter().n];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

Parameter Underlying::get_parameter() const { return PresentationParameter; }

Underlying::Parameter
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

sint16 Underlying::lattice_height() const {
    sint16 n = get_parameter().n;
    return n * (n + 1);
}

Underlying::Underlying(Parameter p)
    : PresentationParameter(p), permutation_table(p.n), CoefficientTable(p.n) {}

void Underlying::print(IndentedOStream &os) const {
    thread_local sint16 dir_perm[MaxN];
    Underlying copy = *this;
    copy.Direct(dir_perm);
    sint16 n = get_parameter().n, e = get_parameter().e;
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

void Underlying::of_string(const std::string &str, size_t &pos) {
    sint16 n = get_parameter().n, e = get_parameter().e;
    std::smatch match;
    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"([st])[\\s\\t]*_?[\\s\\t]*(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[2]);
        pos += match[0].length();
        if (match[1] == "s") {
            if ((3 <= i) && (i <= n)) {
                // s_i.
                identity();
                permutation_table[i - 2] = i - 1;
                permutation_table[i - 1] = i - 2;
            } else {
                throw InvalidStringError(
                    "Invalid index for s type generator!\n" +
                    std::to_string(i) + " is not in [3, " + std::to_string(n) +
                    "].");
            }
        } else {
            // t_i.
            i = Rem(i, e);
            identity();
            permutation_table[0] = 1;
            permutation_table[1] = 0;
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
    for (sint16 i = 0; i < get_parameter().n; i++) {
        dir_perm[permutation_table[i]] = i;
    }
};

Underlying Underlying::left_meet(const Underlying &b) const {
    thread_local sint16 dir_perm_a[MaxN], dir_perm_b[MaxN], dir_perm_meet[MaxN];
    Underlying a_copy = *this;
    Underlying b_copy = b;
    Underlying meet(get_parameter());
    meet.identity();
    a_copy.Direct(dir_perm_a);
    b_copy.Direct(dir_perm_b);
    meet.Direct(dir_perm_meet);
    sint16 n = get_parameter().n, e = get_parameter().e;
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

void Underlying::identity() {
    for (sint16 i = 0; i < get_parameter().n; i++) {
        permutation_table[i] = i;
        CoefficientTable[i] = 0;
    }
}

void Underlying::delta() {
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 1; i < n; i++) {
        permutation_table[i] = i;
        CoefficientTable[i] = 1;
    }
    permutation_table[0] = 0;
    CoefficientTable[0] = Rem(-n + 1, e);
}

bool Underlying::compare(const Underlying &b) const {
    sint16 i;
    for (i = 0; i < get_parameter().n; i++) {
        if ((permutation_table[i] != b.permutation_table[i]) ||
            (CoefficientTable[i] != b.CoefficientTable[i])) {
            return false;
        }
    }
    return true;
};

Underlying Underlying::Inverse() const {
    Underlying f = Underlying(get_parameter());
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 0; i < n; i++) {
        f.permutation_table[permutation_table[i]] = i;
        f.CoefficientTable[permutation_table[i]] =
            CoefficientTable[i] == 0 ? 0 : e - CoefficientTable[i];
    }
    return f;
};

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 0; i < n; i++) {
        f.permutation_table[i] = b.permutation_table[permutation_table[i]];
        f.CoefficientTable[i] = Rem(
            b.CoefficientTable[permutation_table[i]] + CoefficientTable[i], e);
    }
    return f;
};

Underlying Underlying::left_complement(const Underlying &b) const {
    return b.product(Inverse());
};

Underlying Underlying::right_complement(const Underlying &b) const {
    return Inverse().product(b);
};

void Underlying::delta_conjugate_mut(sint16 k) {
    // delta is diagonal, and acts almost homothetically.
    // Therefore conjugating by some power of it does nothing on most
    // coefficients.
    sint16 n = get_parameter().n, e = get_parameter().e;

    // In this case `*this` commutes with delta.
    if (permutation_table[0] == 0) {
        return;
    }

    // Otherwise the two non trivial coefficient are 0 and the i such that
    // `permutation_table[i] == 0`.
    CoefficientTable[0] = Rem(CoefficientTable[0] + k * n, e);
    for (sint16 i = 1; i < n; i++) {
        if (permutation_table[i] == 0) {
            CoefficientTable[i] = Rem(CoefficientTable[i] - k * n, e);
            return;
        }
    }
}

void Underlying::randomize() { throw NonRandomizable(); }

std::vector<Underlying> Underlying::atoms() const {
    Parameter p = get_parameter();
    std::vector<Underlying> atoms;
    Underlying atom = Underlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 2; i <= n - 1; i++) {
        // s_(i+1).
        atom.identity();
        atom.permutation_table[i - 1] = i;
        atom.permutation_table[i] = i - 1;
        atoms.push_back(atom);
    }
    for (sint16 k = 0; k < e; k++) {
        // t_k.
        atom.identity();
        atom.permutation_table[0] = 1;
        atom.permutation_table[1] = 0;
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
    p.print(*this);
    return *this;
};

} // namespace cgarside
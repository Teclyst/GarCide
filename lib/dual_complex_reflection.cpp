#include "dual_complex_reflection.h"
#include <iostream>
#include <string>

namespace cgarside {

namespace dual_complex {

void Parameter::Print(IndentedOStream &os) const {
    os << "(e: " << e << ", n: " << n << ")";
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
        if ((2 <= e) && (2 <= n) && (n <= MaxN) && (e * n <= MaxE * MaxN)) {
            return Parameter(e, n);
        } else if (2 > e) {
            throw InvalidStringError("e should be at least 2!");
            } else if (2 > n) {
            throw InvalidStringError("n should be at least 2!");
        } else if (n > MaxN) {
            throw InvalidStringError("n is too big!\n" +
                                     match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MaxN) + ".");
        } else {
            throw InvalidStringError("e * n is too big!\n" +
                                     match.str(1) + " * " + match.str(2) + " = " + std::to_string(e * n) +
                                     " is strictly greater than " +
                                     std::to_string(MaxN * MaxE) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
}

sint16 Underlying::LatticeHeight() const { return GetParameter().n + 1; }

Underlying::Underlying(Parameter p)
    : PresentationParameter(p), PermutationTable(p.n + 1),
      CoefficientTable(p.n + 1) {}

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

void Underlying::Print(IndentedOStream &os) const {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::vector<bool> seen(n + 1, false);
    std::vector<sint16> curr_cycle;
    sint16 other_smallest = 0, cycle_type, curr, c = 0;
    seen[0] = true;
    bool is_first = true;

    curr = PermutationTable[0];
    // Short assymetric case.
    if (curr != 0) {
        while (curr != 0) {
            curr_cycle.push_back(curr);
            seen[curr] = true;
            other_smallest =
                ((curr < curr_cycle[other_smallest]) && (curr != 0))
                    ? c
                    : other_smallest;
            curr = PermutationTable[curr];
            c++;
        }
        other_smallest =
            ((other_smallest == 0) ? int(curr_cycle.size()) : other_smallest);
        if (int(curr_cycle.size()) > 0) {
            is_first = false;
        }
        for (sint16 i = int(curr_cycle.size()) - 1; i >= 1; i--) {
            os << "("
               << ((i >= other_smallest)
                       ? curr_cycle[i] + Rem(CoefficientTable[0] + 1, e) * n
                       : curr_cycle[i] + CoefficientTable[0] * n)
               << ", "
               << ((i >= other_smallest + 1)
                       ? curr_cycle[i - 1] + Rem(CoefficientTable[0] + 1, e) * n
                       : curr_cycle[i - 1] + CoefficientTable[0] * n)
               << ") ";
        }
        os << curr_cycle[0] + CoefficientTable[0] * n;
    }

    for (sint16 i = 1; i <= n; ++i) {
        if (!seen[i]) {
            seen[i] = true;
            c = 0;
            other_smallest = 0;
            cycle_type = CoefficientTable[i];
            if (CoefficientTable[i] == e - 1) {
                other_smallest = c + 1;
            }
            curr = PermutationTable[i];
            curr_cycle.clear();
            curr_cycle.push_back(i);
            while (curr != i) {
                curr_cycle.push_back(curr);
                seen[curr] = true;
                c++;
                if (CoefficientTable[curr] == e - 1) {
                    other_smallest = c + 1;
                }
                cycle_type += CoefficientTable[curr];
                curr = PermutationTable[curr];
            }
            // If cycle_type == 0, then the cycle is short.
            if (Rem(cycle_type, e) == 0) {
                if (int(curr_cycle.size()) > 1 && !is_first) {
                    os << " ";
                } else if (int(curr_cycle.size()) > 1) {
                    is_first = false;
                }
                for (sint16 i = int(curr_cycle.size()) - 1; i >= 1; i--) {
                    os << "("
                       << ((other_smallest + i >= int(curr_cycle.size()))
                               ? curr_cycle[other_smallest + i -
                                            int(curr_cycle.size())] +
                                     n
                               : curr_cycle[other_smallest + i])
                       << ", "
                       << ((other_smallest + i - 1 >= int(curr_cycle.size()))
                               ? curr_cycle[other_smallest + i - 1 -
                                            int(curr_cycle.size())] +
                                     n
                               : curr_cycle[other_smallest + i - 1])
                       << ")" << ((i == 1) ? "" : " ");
                }
            }
            // Otherwise, it is long.
            else {
                if (!is_first) {
                    os << " ";
                } else {
                    is_first = false;
                }
                other_smallest = ((other_smallest == 0) ? int(curr_cycle.size())
                                                        : other_smallest);
                for (sint16 i = int(curr_cycle.size()) - 1; i >= 1; i--) {
                    os << "("
                       << ((i >= other_smallest) ? curr_cycle[i] + n
                                                 : curr_cycle[i])
                       << ", "
                       << ((i >= other_smallest + 1) ? curr_cycle[i - 1] + n
                                                     : curr_cycle[i - 1])
                       << ") ";
                }
                os << curr_cycle[0] << " " << curr_cycle[0] + n;
            }
        }
    }
}

void Underlying::OfString(const std::string &str, size_t &pos) {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::smatch match;
    if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"D"},
                                 std::regex_constants::match_continuous)) {
        pos += match[0].length();
        Delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"\\([\\s\\t]*(" + number_regex +
                                     ")[\\s\\t]*,?[\\s\\t]*(" + number_regex +
                                     ")[\\s\\t]*\\)"},
                          std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[1]);
        sint16 j = std::stoi(match[2]);
        pos += match[0].length();
        i = Rem(i - 1, e * n) + 1;
        j = Rem(j - 1, e * n) + 1;
        if (i > j) {
            std::swap(i, j);
        }
        if ((i < n) && (j > (e - 1) * n)) {
            std::swap(i, j);
            j += e * n;
        }

        if ((j - i < n) && (j != i)) {
            i = Rem(i - 1, n) + 1;
            j = Rem(j - 1, n) + 1;
            j = (j < i) ? j + n : j;
            Identity();
            PermutationTable[i] = (j > n) ? j - n : j;
            PermutationTable[(j > n) ? j - n : j] = i;
            CoefficientTable[i] = (j > n) ? 1 : 0;
            CoefficientTable[(j > n) ? j - n : j] = (j > n) ? e - 1 : 0;
        } else if ((Rem(i, n) == Rem(j, n))) {
            throw InvalidStringError("Indexes for short symmetric generators "
                                     "should not be equal mod " +
                                     std::to_string(e * n) + "!\n(" +
                                     std::to_string(std::stoi(match[1])) +
                                     ", " +
                                     std::to_string(std::stoi(match[2])) +
                                     ") is not a valid factor.");
        } else {
            throw InvalidStringError(
                "Indexes for short generators should be at most " +
                std::to_string(n - 1) + " apart mod" + std::to_string(e * n) +
                "!\n(" + std::to_string(std::stoi(match[1])) + ", " +
                std::to_string(std::stoi(match[2])) +
                ") is not a valid factor.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"(" + number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[1]);
        pos += match[0].length();
        i = Rem(i - 1, e * n) + 1;
        sint16 q = Rem(Quot(i - 1, n), e);
        sint16 r = Rem(i - 1, n) + 1;
        Identity();
        PermutationTable[0] = r;
        PermutationTable[r] = 0;
        CoefficientTable[0] = q;
        CoefficientTable[r] = (q == 0) ? 0 : e - q;
    } else {
        throw InvalidStringError(
            "Could not extract a factor from \"" + str.substr(pos) +
            "\"!\nA factor should match regex '(' Z ','? Z ')' | Z | "
            "'D',\nwhere Z matches integers, and ignoring whitespaces.");
    }
}

void Underlying::AssignPartition(sint16 *x) const {
    sint16 n = GetParameter().n, e = GetParameter().e;
    std::vector<sint16> curr_cycle;
    sint16 other_smallest, cycle_type, curr, c;
    for (sint16 i = 1; i <= e * n; ++i) {
        x[i] = -1;
    }
    x[0] = 0;

    curr = PermutationTable[0];
    if (curr != 0) {
        other_smallest = 0;
        c = 0;
        while (curr != 0) {
            curr_cycle.push_back(curr);
            other_smallest =
                ((curr < curr_cycle[other_smallest]) && (curr != 0))
                    ? c
                    : other_smallest;
            curr = PermutationTable[curr];
            c++;
        }
        if (other_smallest != 0) {
            for (sint16 k = 0; k < other_smallest; k++) {
                for (sint16 l = 0; l < e - 1; l++) {
                    x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
                }
                x[curr_cycle[k] + (e - 1) * n] = curr_cycle[other_smallest];
                x[curr_cycle[k] + Rem(CoefficientTable[0], e) * n] = 0;
            }
            for (sint16 k = other_smallest; k < int(curr_cycle.size()); k++) {
                for (sint16 l = 0; l < e - 1; l++) {
                    x[curr_cycle[k] + (l + 1) * n] = curr_cycle[0] + l * n;
                }
                x[curr_cycle[k]] = curr_cycle[other_smallest];
                x[curr_cycle[k] + Rem(CoefficientTable[0] + 1, e) * n] = 0;
            }
        } else {
            for (sint16 k = 0; k < int(curr_cycle.size()); k++) {
                for (sint16 l = 0; l < e; l++) {
                    x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
                }
                x[curr_cycle[k] + Rem(CoefficientTable[0], e) * n] = 0;
            }
        }
    }

    for (sint16 i = 1; i <= n; ++i) {
        if (x[i] < 0) {
            c = 0;
            other_smallest = 0;
            cycle_type = CoefficientTable[i];
            if (CoefficientTable[i] == e - 1) {
                other_smallest = 1;
            }
            curr = PermutationTable[i];
            curr_cycle.clear();
            curr_cycle.push_back(i);
            while (curr != i) {
                curr_cycle.push_back(curr);
                c++;
                if (CoefficientTable[curr] == e - 1) {
                    other_smallest = c + 1;
                }
                cycle_type += CoefficientTable[curr];
                curr = PermutationTable[curr];
            }
            // If cycle_type == 0, then the cycle is short.
            if (Rem(cycle_type, e) == 0) {
                if (other_smallest != 0) {
                    for (sint16 k = 0; k < other_smallest; k++) {
                        x[curr_cycle[k]] = i;
                        for (sint16 l = 0; l < e - 1; l++) {
                            x[curr_cycle[k] + (l + 1) * n] =
                                curr_cycle[other_smallest] + l * n;
                        }
                    }
                    for (sint16 k = other_smallest; k < int(curr_cycle.size());
                         k++) {
                        x[curr_cycle[k] + (e - 1) * n] = i;
                        for (sint16 l = 0; l < e - 1; l++) {
                            x[curr_cycle[k] + l * n] =
                                curr_cycle[other_smallest] + l * n;
                        }
                    }
                } else {
                    for (sint16 k = 0; k < int(curr_cycle.size()); k++) {
                        for (sint16 l = 0; l < e; l++) {
                            x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
                        }
                    }
                }
            }
            // Otherwise, it is long.
            else {
                for (sint16 k = 0; k < int(curr_cycle.size()); k++) {
                    for (sint16 l = 0; l < e; l++) {
                        x[curr_cycle[k] + l * n] = 0;
                    }
                }
            }
        }
    }
}

// We assume e > 1.
void Underlying::OfPartition(const sint16 *x) {
    thread_local sint16 z[MaxN + 1];
    sint16 min_cycle_0 = 0, max_cycle_0 = 0;
    sint16 n = GetParameter().n, e = GetParameter().e, r;

    for (sint16 i = 0; i <= n; ++i) {
        z[i] = -1;
        PermutationTable[i] = -1;
        CoefficientTable[i] = -1;
    }
    // First find short symmetrical cycles.
    for (sint16 i = 2 * n - 1; i >= n + 1; --i) {
        if ((x[i] <= n) && (x[i] >= 1)) {
            r = i - n;
            if (z[x[i]] == -1) {
                PermutationTable[r] = x[i];
                CoefficientTable[r] = e - 1;
                z[x[i]] = r;
            } else {
                PermutationTable[r] = z[x[i]];
                CoefficientTable[r] = 0;
                z[x[i]] = r;
            }
        }
    }
    for (sint16 i = n; i >= 1; --i) {
        if ((x[i] <= n) && (x[i] >= 1) && (x[i + n] > n)) {
            if ((z[x[i]] == -1)) {
                PermutationTable[i] = x[i];
                CoefficientTable[i] = 0;
            } else {
                CoefficientTable[i] = (z[x[i]] < i) ? 1 : 0;
                PermutationTable[i] = z[x[i]];
            }
            z[x[i]] = i;
        }
    }

    // Then consider the part with 0 in it.
    for (sint16 i = 1; i <= e * n; i++) {
        if (x[i] == 0) {
            min_cycle_0 = i;
            break;
        }
    }
    if (min_cycle_0 != 0) {
        // Determine if it is long symmetric.
        if (x[Rem(min_cycle_0 + n - 1, e * n) + 1] == 0) {
            CoefficientTable[0] = e - 1;
            PermutationTable[0] = 0;
            for (sint16 i = n; i >= 1; i--) {
                if (x[i] == 0) {
                    if (z[x[i]] == -1) {
                        PermutationTable[i] = min_cycle_0;
                        CoefficientTable[i] = 1;
                    } else {
                        PermutationTable[i] = z[x[i]];
                        CoefficientTable[i] = 0;
                    }
                    z[x[i]] = i;
                }
            }
        }
        // Assymetric case.
        else {
            for (sint16 i = e * n; i >= 1; i--) {
                if (x[i] == 0) {
                    max_cycle_0 = i;
                    break;
                }
            }
            if ((min_cycle_0 <= n) && (max_cycle_0 > (e - 1) * n)) {
                for (sint16 i = (e - 1) * n + 1; i <= e * n; i++) {
                    if (x[i] == 0) {
                        min_cycle_0 = i;
                        break;
                    }
                }
                for (sint16 i = n; i >= 1; i--) {
                    if (x[i] == 0) {
                        max_cycle_0 = i;
                        break;
                    }
                }
            }
            sint16 q_min = Quot(min_cycle_0 - 1, n),
                   q_max = Quot(max_cycle_0 - 1, n);
            sint16 r_min = Rem(min_cycle_0 - 1, n) + 1;
            z[0] = 0;
            CoefficientTable[0] = q_min;
            PermutationTable[0] = r_min;

            for (sint16 i = n - 1; i >= n - r_min + 1; --i) {
                sint16 i_en = Rem(i + min_cycle_0 - 1, e * n) + 1;
                r = Rem(i + min_cycle_0 - 1, n) + 1;
                if (x[i_en] == 0) {
                    PermutationTable[r] = z[0];
                    CoefficientTable[r] = (z[0] == 0) ? Rem(e - q_max, e) : 0;
                    z[0] = r;
                }
            }
            for (sint16 i = n - r_min; i >= 0; --i) {
                sint16 i_en = Rem(i + min_cycle_0 - 1, e * n) + 1;
                r = Rem(i + min_cycle_0 - 1, n) + 1;
                if ((x[i_en] == 0)) {
                    if (z[0] == 0) {
                        PermutationTable[r] = z[0];
                        CoefficientTable[r] = Rem(e - q_min, e);
                    } else {
                        CoefficientTable[r] = (z[0] < r) ? 1 : 0;
                        PermutationTable[r] = z[0];
                    }
                    z[0] = r;
                }
            }
        }
    }

    // Short symmetric case.
    else {
        CoefficientTable[0] = 0;
        PermutationTable[0] = 0;
    }
}

Underlying Underlying::LeftMeet(const Underlying &b) const {
    thread_local sint16 x[MaxE * MaxN + 1], y[MaxE * MaxN + 1],
        z[MaxE * MaxN + 1];

    AssignPartition(x);
    b.AssignPartition(y);

    thread_local sint16 P[MaxE * MaxN + 1][MaxE * MaxN + 1];

    for (sint16 i = GetParameter().e * GetParameter().n; i >= 0; i--) {
        P[x[i]][y[i]] = i;
    }

    for (sint16 i = 0; i <= GetParameter().e * GetParameter().n; i++) {
        z[i] = P[x[i]][y[i]];
    }

    Underlying c = Underlying(*this);

    c.OfPartition(z);

    return c;
}

void Underlying::Identity() {
    for (sint16 i = 0; i <= GetParameter().n; i++) {
        PermutationTable[i] = i;
        CoefficientTable[i] = 0;
    }
}

void Underlying::Delta() {
    sint16 i, n = GetParameter().n;
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = i + 1;
        CoefficientTable[i] = 0;
    }
    PermutationTable[0] = 0;
    CoefficientTable[0] = GetParameter().e - 1;
    PermutationTable[n] = 1;
    CoefficientTable[n] = 1;
}

bool Underlying::Compare(const Underlying &b) const {
    sint16 i;
    for (i = 0; i <= GetParameter().n; i++) {
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
    for (i = 0; i <= n; i++) {
        if ((PermutationTable[i] > n) || (PermutationTable[i] < 0)) {
            while (true) {
            };
        }
        f.PermutationTable[PermutationTable[i]] = i;
        f.CoefficientTable[PermutationTable[i]] =
            CoefficientTable[i] == 0 ? 0 : e - CoefficientTable[i];
    }
    return f;
};

Underlying Underlying::Product(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    for (i = 0; i <= n; i++) {
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
    sint16 i, n = GetParameter().n, e = GetParameter().e;
    sint16 q = Quot(k, n), r = Rem(k, n), q_e = Rem(q, e);

    Underlying delta_k = Underlying(GetParameter());

    delta_k.PermutationTable[0] = 0;
    delta_k.CoefficientTable[0] = Rem(-k, e);
    for (i = 1; i <= n - r; i++) {
        delta_k.PermutationTable[i] = Rem(i + k - 1, n) + 1;
        delta_k.CoefficientTable[i] = q_e;
    }
    q_e += 1;
    q_e = ((q_e == e) ? 0 : q_e);
    for (i = n - r + 1; i <= n; i++) {
        delta_k.PermutationTable[i] = Rem(i + k - 1, n) + 1;
        delta_k.CoefficientTable[i] = q_e;
    }

    *this = delta_k.Inverse().Product((*this).Product(delta_k));
}

void Underlying::Randomize() { throw NonRandomizable(); }

std::vector<Underlying> Underlying::Atoms() const {
    Parameter p = GetParameter();
    std::vector<Underlying> atoms;
    Underlying atom = Underlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 1; i <= n; i++) {
        for (sint16 j = i + 1; j <= n; j++) {
            atom.Identity();
            atom.PermutationTable[i] = j;
            atom.PermutationTable[j] = i;
            atoms.push_back(atom);
        }
        for (sint16 j = 1; j < i; j++) {
            atom.Identity();
            atom.PermutationTable[i] = j;
            atom.PermutationTable[j] = i;
            atom.CoefficientTable[i] = 1;
            atom.CoefficientTable[j] = e - 1;
            atoms.push_back(atom);
        }
    }
    for (sint16 i = 1; i <= n; i++) {
        atom.Identity();
        atom.PermutationTable[0] = i;
        atom.PermutationTable[i] = 0;
        atom.CoefficientTable[0] = 0;
        atom.CoefficientTable[i] = 0;
        atoms.push_back(atom);
    }
    for (sint16 k = 1; k < e; k++) {
        for (sint16 i = 1; i <= n; i++) {
            atom.Identity();
            atom.PermutationTable[0] = i;
            atom.PermutationTable[i] = 0;
            atom.CoefficientTable[0] = k;
            atom.CoefficientTable[i] = e - k;
            atoms.push_back(atom);
        }
    }
    return atoms;
}

std::size_t Underlying::Hash() const {
    std::size_t h = 0;
    for (sint16 i = 1; i <= GetParameter().n; i++) {
        h = h * 31 + PermutationTable[i];
    }
    for (sint16 i = 1; i <= GetParameter().n; i++) {
        h = h * 31 + CoefficientTable[i];
    }
    return h;
}

} // namespace dual_complex

template <>
IndentedOStream &
    IndentedOStream::operator<< <dual_complex::Parameter>(const dual_complex::Parameter &p) {
    p.Print(*this);
    return *this;
};

} // namespace cgarside
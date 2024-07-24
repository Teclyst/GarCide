#include "band_braid.h"
#include "cbraid.h"
#include <iostream>
#include <string>

namespace cgarside::band {

sint16 Underlying::GetParameter() const { return PresentationParameter; }

sint16 Underlying::LatticeHeight() const { return GetParameter(); }

Underlying::ParameterType
Underlying::parameter_of_string(const std::string &str) {
    std::smatch match;

    if (std::regex_match(str, match,
                         std::regex{"[\\s\\t]*(" + number_regex + ")[\\s\\t]*"},
                         std::regex_constants::match_continuous)) {
        sint16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Number of strands is too big!\n" +
                                     match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        if (((2 <= i) && (i <= MaxBraidIndex))) {
            return i;
        } else if (2 > i) {
            throw InvalidStringError("Number of strands should be at least 2!");
        } else {
            throw InvalidStringError("Number of strands is too big!\n" +
                                     match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MaxBraidIndex) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

Underlying::Underlying(sint16 n)
    : PresentationParameter(n), PermutationTable(n + 1) {}

void Underlying::Print(IndentedOStream &os) const {
    // Recall that a band braid is represented by decreasing cycles.
    sint16 i, j, n = GetParameter();
    std::vector<sint16> curr_cycle;
    std::vector<bool> seen(n + 1, false);
    bool is_first = true;
    for (i = 1; i <= n; ++i) {
        if (!seen[i]) {
            j = i;
            curr_cycle.clear();
            while (j < PermutationTable[j]) {
                curr_cycle.push_back(j);
                seen[j] = true;
                j = PermutationTable[j];
            }
            curr_cycle.push_back(j);
            seen[j] = true;
            if (int(curr_cycle.size()) > 1 && !is_first) {
                os << " ";
            } else if (int(curr_cycle.size()) > 1) {
                is_first = false;
            }
            for (j = int(curr_cycle.size()) - 1; j >= 1; --j) {
                os << "(" << curr_cycle[j] << ", " << curr_cycle[j - 1] << ")"
                   << ((j == 1) ? "" : " ");
            }
        }
    }
}

void Underlying::OfString(const std::string &str, size_t &pos) {
    sint16 n = GetParameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        Delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"\\([\\s\\t]*(" + number_regex +
                                            ")[\\s\\t]*,?[\\s\\t]*(" +
                                            number_regex + ")[\\s\\t]*\\)"},
                                 std::regex_constants::match_continuous)) {
        pos += match[0].length();
        sint16 i, j;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        try {
            j = std::stoi(match[2]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(2) +
                                     " can not be converted to a C++ integer.");
        }
        if ((i >= 1) && (i <= n) && (j >= 1) && (j <= n) && (i != j)) {
            Identity();
            PermutationTable[i] = j;
            PermutationTable[j] = i;
        } else if ((i < 1) || (i > n)) {
            throw InvalidStringError("Invalid index for dual generator!\n" +
                                     match.str(1) + " is not in [1, " +
                                     std::to_string(n) + "].");
        } else if ((j < 1) || (j > n)) {
            throw InvalidStringError("Invalid index for dual generator!\n" +
                                     match.str(2) + " is not in [1, " +
                                     std::to_string(n) + "].");
        } else {
            throw InvalidStringError(
                "Indexes for dual generators should not be equal!\n(" +
                match.str(1) + ", " + match.str(2) +
                ") is not a valid factor.");
        }
    } else {
        throw InvalidStringError(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex '('Z ','? Z ')' | 'D',\nwhere Z "
            "matches integers, and ignoring whitespaces.");
    }
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
    for (sint16 i = 1; i < GetParameter(); i++) {
        os << PermutationTable[i] << ", ";
    }
    os << PermutationTable[GetParameter()];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

void Underlying::AssignDCDT(sint16 *x) const {
    for (sint16 i = 1; i <= GetParameter(); ++i)
        x[i] = 0;
    for (sint16 i = 1; i <= GetParameter(); ++i) {
        if (x[i] == 0)
            x[i] = i;
        if (PermutationTable[i] > i)
            x[PermutationTable[i]] = x[i];
    }
}

void Underlying::OfDCDT(const sint16 *x) {
    thread_local sint16 z[MaxBraidIndex];

    for (sint16 i = 1; i <= GetParameter(); ++i)
        z[i] = 0;
    for (sint16 i = GetParameter(); i >= 1; --i) {
        PermutationTable[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
        z[x[i]] = i;
    }
}

Underlying Underlying::LeftMeet(const Underlying &b) const {
    thread_local sint16 x[MaxBraidIndex], y[MaxBraidIndex], z[MaxBraidIndex];

    AssignDCDT(x);
    b.AssignDCDT(y);

    thread_local sint16 P[MaxBraidIndex][MaxBraidIndex];

    for (sint16 i = GetParameter(); i >= 1; i--) {
        P[x[i]][y[i]] = i;
    }

    for (sint16 i = 1; i <= GetParameter(); i++) {
        z[i] = P[x[i]][y[i]];
    }

    Underlying c = Underlying(*this);

    c.OfDCDT(z);

    return c;
}

Underlying Underlying::RightMeet(const Underlying &b) const {
    return LeftMeet(b);
}

void Underlying::Identity() {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = i;
    }
}

void Underlying::Delta() {
    sint16 i, n = GetParameter();
    for (i = 1; i < n; i++) {
        PermutationTable[i] = i + 1;
    }
    PermutationTable[n] = 1;
}

bool Underlying::Compare(const Underlying &b) const {
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        if (PermutationTable[i] != b.PermutationTable[i]) {
            return false;
        }
    }
    return true;
};

Underlying Underlying::Inverse() const {
    Underlying f = Underlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
};

Underlying Underlying::Product(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
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
    Underlying under = *this;
    sint16 i, n = GetParameter();

    for (i = 1; i <= n; i++) {
        under.PermutationTable[i] =
            Rem(PermutationTable[Rem(i - k - 1, n) + 1] + k - 1, n) + 1;
    }
    *this = under;
}

void Underlying::Randomize() {
    CBraid::BandPresentation pres = CBraid::BandPresentation(GetParameter());
    pres.Randomize(PermutationTable.data());
}

std::vector<Underlying> Underlying::Atoms() const {
    sint16 n = GetParameter();
    std::vector<Underlying> atoms;
    for (sint16 i = 1; i <= n; i++) {
        for (sint16 j = 1; j < i; j++) {
            Underlying atom = Underlying(n);
            atom.Identity();
            atom.PermutationTable[i] = j;
            atom.PermutationTable[j] = i;
            atoms.push_back(atom);
        }
    }
    return atoms;
}

size_t Underlying::Hash() const {
    std::size_t h = 0;
    for (sint16 i = 1; i <= GetParameter(); i++) {
        h = h * 31 + PermutationTable[i];
    }
    return h;
}

} // namespace cgarside::band

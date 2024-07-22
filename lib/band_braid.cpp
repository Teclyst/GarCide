#include "band_braid.h"
#include "cbraid.h"
#include <iostream>
#include <string>

namespace cgarside {

sint16 BandBraidUnderlying::GetParameter() const {
    return PresentationParameter;
}

sint16 BandBraidUnderlying::LatticeHeight() const { return GetParameter(); }

BandBraidUnderlying::BandBraidUnderlying(sint16 n)
    : PresentationParameter(n), PermutationTable(n + 1) {}

void BandBraidUnderlying::Print(IndentedOStream &os) const {
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

void BandBraidUnderlying::OfString(const std::string &str, size_t &pos) {
    sint16 n = GetParameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"\\([\\s\\t]*(" + number_regex +
                                     ")[\\s\\t]*,?[\\s\\t]*(" + number_regex +
                                     ")[\\s\\t]*\\)"},
                          std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[1]);
        sint16 j = std::stoi(match[2]);
        pos += match[0].length();
        if ((i >= 1) && (i <= n) && (j >= 1) && (j <= n) && (i != j)) {
            Identity();
            PermutationTable[i] = j;
            PermutationTable[j] = i;
        } else if ((i < 1) || (i > n)) {
            throw InvalidStringError("Invalid index for dual generator!\n" +
                                     std::to_string(i) + " is not in [1, " +
                                     std::to_string(n) + "].");
        } else if ((j < 1) || (j > n)) {
            throw InvalidStringError("Invalid index for dual generator!\n" +
                                     std::to_string(j) + " is not in [1, " +
                                     std::to_string(n) + "].");
        } else {
            throw InvalidStringError(
                "Indexes for dual generators should not be equal!\n(" +
                std::to_string(i) + ", " + std::to_string(j) +
                ") is not a valid factor.");
        }
    } else {
        throw InvalidStringError(
            "Could not extract a factor from \"" + str.substr(pos) +
            "\"!\nA factor should match regex \\([1 - 9] [0 - 9]*,? [1 - 9] [0 "
            "- 9]*\\) (ignoring whitespaces).");
    }
}

void BandBraidUnderlying::Debug(IndentedOStream &os) const {
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

void BandBraidUnderlying::AssignDCDT(sint16 *x) const {
    for (sint16 i = 1; i <= GetParameter(); ++i)
        x[i] = 0;
    for (sint16 i = 1; i <= GetParameter(); ++i) {
        if (x[i] == 0)
            x[i] = i;
        if (PermutationTable[i] > i)
            x[PermutationTable[i]] = x[i];
    }
}

void BandBraidUnderlying::OfDCDT(const sint16 *x) {
    thread_local sint16 z[MaxBraidIndex];

    for (sint16 i = 1; i <= GetParameter(); ++i)
        z[i] = 0;
    for (sint16 i = GetParameter(); i >= 1; --i) {
        PermutationTable[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
        z[x[i]] = i;
    }
}

BandBraidUnderlying
BandBraidUnderlying::LeftMeet(const BandBraidUnderlying &b) const {
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

    BandBraidUnderlying c = BandBraidUnderlying(*this);

    c.OfDCDT(z);

    return c;
}

BandBraidUnderlying
BandBraidUnderlying::RightMeet(const BandBraidUnderlying &b) const {
    return LeftMeet(b);
}

void BandBraidUnderlying::Identity() {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = i;
    }
}

void BandBraidUnderlying::Delta() {
    sint16 i, n = GetParameter();
    for (i = 1; i < n; i++) {
        PermutationTable[i] = i + 1;
    }
    PermutationTable[n] = 1;
}

bool BandBraidUnderlying::Compare(const BandBraidUnderlying &b) const {
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        if (PermutationTable[i] != b.PermutationTable[i]) {
            return false;
        }
    }
    return true;
};

BandBraidUnderlying BandBraidUnderlying::Inverse() const {
    BandBraidUnderlying f = BandBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
};

BandBraidUnderlying
BandBraidUnderlying::Product(const BandBraidUnderlying &b) const {
    BandBraidUnderlying f = BandBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
    }
    return f;
};

BandBraidUnderlying
BandBraidUnderlying::LeftComplement(const BandBraidUnderlying &b) const {
    return b.Product(Inverse());
};

BandBraidUnderlying
BandBraidUnderlying::RightComplement(const BandBraidUnderlying &b) const {
    return Inverse().Product(b);
};

void BandBraidUnderlying::DeltaConjugate(sint16 k) {
    BandBraidUnderlying under = *this;
    sint16 i, n = GetParameter();

    for (i = 1; i <= n; i++) {
        under.PermutationTable[i] =
            Rem(PermutationTable[Rem(i - k - 1, n) + 1] + k - 1, n) + 1;
    }
    *this = under;
}

void BandBraidUnderlying::Randomize() {
    CBraid::BandPresentation pres = CBraid::BandPresentation(GetParameter());
    pres.Randomize(PermutationTable.data());
}

std::vector<BandBraidUnderlying> BandBraidUnderlying::Atoms() const {
    sint16 n = GetParameter();
    std::vector<BandBraidUnderlying> atoms;
    for (sint16 i = 1; i <= n; i++) {
        for (sint16 j = 1; j < i; j++) {
            BandBraidUnderlying atom = BandBraidUnderlying(n);
            atom.Identity();
            atom.PermutationTable[i] = j;
            atom.PermutationTable[j] = i;
            atoms.push_back(atom);
        }
    }
    return atoms;
}

size_t BandBraidUnderlying::Hash() const {
        std::size_t h = 0;
        for (sint16 i = 1; i <= GetParameter(); i++) {
            h = h * 31 + PermutationTable[i];
        }
        return h;
    }

} // namespace CGarside

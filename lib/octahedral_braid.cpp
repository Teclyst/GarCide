#include "octahedral_braid.h"
#include <iostream>
#include <string>

namespace cgarside::octahedral {

/**
 * @brief Maximum braid index.
 *
 * The greatest index that may be used for braids.
 *
 * It is used because we use `thread_local` objects to avoid some allocations,
 * and their size must be known at compile time.
 *
 * Having too big `thread_local` objects might cause some issue with thread
 * spawning.
 */
const sint16 MaxBraidIndex = 256;

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
            throw InvalidStringError("Parameter is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        if (((1 <= i) && (i <= MaxBraidIndex))) {
            return i;
        } else if (1 > i) {
            throw InvalidStringError("Parameter should be at least 1!");
        } else {
            throw InvalidStringError("Parameter strands is too big!\n" +
                                     match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MaxBraidIndex) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

sint16 Underlying::GetParameter() const { return PresentationParameter; }

sint16 Underlying::LatticeHeight() const { return GetParameter(); }

Underlying::Underlying(sint16 n)
    : PresentationParameter(n), PermutationTable(2 * n + 1) {}

void Underlying::Print(IndentedOStream &os) const {
    // Recall that a band braid is represented by decreasing cycles.
    sint16 i, j, n = GetParameter();
    std::vector<sint16> curr_cycle;
    std::vector<bool> seen(n + 1, false);
    for (i = 1; i <= n; i++) {
        seen[i] = false;
    }

    bool is_first = true;

    for (i = 1; i <= n; ++i) {
        if (!seen[i]) {
            curr_cycle.clear();
            j = i;
            while (!seen[Rem(j - 1, n) + 1]) {
                curr_cycle.push_back(j);
                seen[Rem(j - 1, n) + 1] = true;
                j = PermutationTable[j];
            }
            if (!is_first && int(curr_cycle.size()) > 1) {
                os << " ";
            } else if (int(curr_cycle.size()) > 1) {
                is_first = false;
            }
            for (sint16 l = int(curr_cycle.size()) - 1; l >= 1; --l) {
                os << "(" << curr_cycle[l] << ", " << curr_cycle[l - 1] << ")"
                   << ((l == 1) ? "" : " ");
            }
            // Long cycle.
            if (j > n) {
                os << (is_first ? "" : " ") << j - n;
            }
        }
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
    for (sint16 i = 1; i < 2 * GetParameter(); i++) {
        os << PermutationTable[i] << ", ";
    }
    os << PermutationTable[2 * GetParameter()];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
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
        i = Rem(i - 1, 2 * n) + 1;
        j = Rem(j - 1, 2 * n) + 1;
        pos += match[0].length();
        if ((i != j) && (i != Rem(j + n - 1, 2 * n) + 1)) {
            Identity();
            PermutationTable[i] = j;
            PermutationTable[j] = i;
            PermutationTable[Rem(i + n - 1, 2 * n) + 1] =
                Rem(j + n - 1, 2 * n) + 1;
            PermutationTable[Rem(j + n - 1, 2 * n) + 1] =
                Rem(i + n - 1, 2 * n) + 1;
        } else {
            throw InvalidStringError(
                "Indexes for short generators should not be equal mod " +
                std::to_string(n) + "!\n(" + match.str(1) + ", " +
                match.str(2) + ") is not a valid factor.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"(" + number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        sint16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        i = Rem(i - 1, 2 * n) + 1;
        pos += match[0].length();
        Identity();
        PermutationTable[i] = Rem(i + n - 1, 2 * n) + 1;
        PermutationTable[Rem(i + n - 1, 2 * n) + 1] = i;
    } else {
        throw InvalidStringError(
            "Could not extract a factor from \"" + str.substr(pos) +
            "\"!\nA factor should match regex '(' Z ','? Z ')' | Z | "
            "'D',\nwhere Z matches integers, and ignoring whitespaces.");
    }
}

void Underlying::AssignDCDT(sint16 *x) const {
    for (sint16 i = 1; i <= 2 * GetParameter(); ++i)
        x[i] = 0;
    for (sint16 i = 1; i <= 2 * GetParameter(); ++i) {
        if (x[i] == 0)
            x[i] = i;
        if (PermutationTable[i] > i)
            x[PermutationTable[i]] = x[i];
    }
}

void Underlying::OfDCDT(const sint16 *x) {
    thread_local sint16 z[2 * MaxBraidIndex + 1];

    for (sint16 i = 1; i <= 2 * GetParameter(); ++i)
        z[i] = 0;
    for (sint16 i = 2 * GetParameter(); i >= 1; --i) {
        PermutationTable[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
        z[x[i]] = i;
    }
}

Underlying Underlying::LeftMeet(const Underlying &b) const {
    thread_local sint16 x[2 * MaxBraidIndex + 1], y[2 * MaxBraidIndex + 1],
        z[2 * MaxBraidIndex + 1];

    AssignDCDT(x);
    b.AssignDCDT(y);

    thread_local sint16 P[2 * MaxBraidIndex + 1][2 * MaxBraidIndex + 1];

    for (sint16 i = 2 * GetParameter(); i >= 1; i--) {
        P[x[i]][y[i]] = i;
    }

    for (sint16 i = 1; i <= 2 * GetParameter(); i++) {
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
    for (sint16 i = 1; i <= 2 * GetParameter(); i++) {
        PermutationTable[i] = i;
    }
}

void Underlying::Delta() {
    sint16 i, n = GetParameter();
    for (i = 1; i < 2 * n; i++) {
        PermutationTable[i] = i + 1;
    }
    PermutationTable[2 * n] = 1;
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
    for (i = 1; i <= 2 * GetParameter(); i++) {
        f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
};

Underlying Underlying::Product(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    sint16 i;
    for (i = 1; i <= 2 * GetParameter(); i++) {
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

    for (i = 1; i <= 2 * n; i++) {
        under.PermutationTable[i] =
            Rem(PermutationTable[Rem(i - k - 1, 2 * n) + 1] + k - 1, 2 * n) + 1;
    }
    *this = under;
}

void Underlying::Randomize() { throw NonRandomizable(); }

std::vector<Underlying> Underlying::Atoms() const {
    sint16 n = GetParameter();
    std::vector<Underlying> atoms;
    for (sint16 i = 1; i <= n; i++) {
        for (sint16 j = i + 1; j <= n; j++) {
            Underlying atom = Underlying(n);
            atom.Identity();
            atom.PermutationTable[i] = j;
            atom.PermutationTable[j] = i;
            atom.PermutationTable[i + n] = j + n;
            atom.PermutationTable[j + n] = i + n;
            atoms.push_back(atom);
        }
        for (sint16 j = n + i + 1; j <= 2 * n; j++) {
            Underlying atom = Underlying(n);
            atom.Identity();
            atom.PermutationTable[i] = j;
            atom.PermutationTable[j] = i;
            atom.PermutationTable[i + n] = j - n;
            atom.PermutationTable[j - n] = i + n;
            atoms.push_back(atom);
        }
        Underlying atom = Underlying(n);
        atom.Identity();
        atom.PermutationTable[i] = n + i;
        atom.PermutationTable[n + i] = i;
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace cgarside::octahedral
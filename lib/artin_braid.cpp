#include "artin_braid.h"

namespace CGarside {

void ArtinBraidUnderlying::MeetSub(const sint16 *a, const sint16 *b, sint16 *r,
                                   sint16 s, sint16 t) {
    thread_local sint16 u[MaxBraidIndex], v[MaxBraidIndex], w[MaxBraidIndex];

    if (s >= t)
        return;
    sint16 m = (s + t) / 2;
    MeetSub(a, b, r, s, m);
    MeetSub(a, b, r, m + 1, t);

    u[m] = a[r[m]];
    v[m] = b[r[m]];
    if (s < m) {
        for (sint16 i = m - 1; i >= s; --i) {
            u[i] = std::min(a[r[i]], u[i + 1]);
            v[i] = std::min(b[r[i]], v[i + 1]);
        }
    }
    u[m + 1] = a[r[m + 1]];
    v[m + 1] = b[r[m + 1]];
    if (t > m + 1) {
        for (sint16 i = m + 2; i <= t; ++i) {
            u[i] = std::max(a[r[i]], u[i - 1]);
            v[i] = std::max(b[r[i]], v[i - 1]);
        }
    }

    sint16 p = s;
    sint16 q = m + 1;
    for (sint16 i = s; i <= t; ++i)
        w[i] = ((p > m) || (q <= t && u[p] > u[q] && v[p] > v[q])) ? r[q++]
                                                                   : r[p++];
    for (sint16 i = s; i <= t; ++i)
        r[i] = w[i];
}

sint16 ArtinBraidUnderlying::GetParameter() const {
    return PresentationParameter;
};

ArtinBraidUnderlying::ArtinBraidUnderlying(sint16 n)
    : PresentationParameter(n), PermutationTable(n + 1) {}

void ArtinBraidUnderlying::OfString(const std::string &str, size_t &pos) {
    sint16 n = GetParameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"(" + number_regex + ")"},
                          std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[1]);
        pos += match[0].length();
        if ((i >= 1) && (i < n)) {
            Identity();
            PermutationTable[i] = i + 1;
            PermutationTable[i + 1] = i;
        } else {
            throw InvalidStringError("Invalid index for Artin generator!\n" +
                                     std::to_string(i) + " is not in [1, " +
                                     std::to_string(n) + "[.");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract a factor from \"" + str.substr(pos) +
                        "\"!\nA factor should match regex [1 - 9] [0 - 9]*."));
    }
};

void ArtinBraidUnderlying::Print(IndentedOStream &os) const {
    sint16 i, j, k, n = GetParameter();

    ArtinBraidUnderlying c = ArtinBraidUnderlying(*this);

    bool is_first = true;

    for (i = 2; i <= n; i++) {
        for (j = i; j > 1 && c.PermutationTable[j] < c.PermutationTable[j - 1];
             j--) {
            os << (is_first ? "" : " ") << j - 1;
            is_first = false;
            k = c.PermutationTable[j];
            c.PermutationTable[j] = c.PermutationTable[j - 1];
            c.PermutationTable[j - 1] = k;
        }
    }
};

void ArtinBraidUnderlying::Debug(IndentedOStream &os) const {
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

sint16 ArtinBraidUnderlying::LatticeHeight() const {
    sint16 n = GetParameter();
    return n * (n - 1) / 2;
}

void ArtinBraidUnderlying::Identity() {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = i;
    }
};

void ArtinBraidUnderlying::Delta() {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = n + 1 - i;
    }
};

ArtinBraidUnderlying
ArtinBraidUnderlying::LeftMeet(const ArtinBraidUnderlying &b) const {
    thread_local sint16 s[MaxBraidIndex];

    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());

    for (sint16 i = 1; i <= GetParameter(); ++i)
        s[i] = i;
    MeetSub(PermutationTable.data(), b.PermutationTable.data(), s, 1,
            GetParameter());
    for (sint16 i = 1; i <= GetParameter(); ++i)
        f.PermutationTable[s[i]] = i;

    return f;
};

ArtinBraidUnderlying
ArtinBraidUnderlying::RightMeet(const ArtinBraidUnderlying &b) const {
    thread_local sint16 u[MaxBraidIndex], v[MaxBraidIndex];

    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());

    for (sint16 i = 1; i <= GetParameter(); ++i) {
        u[PermutationTable[i]] = i;
        v[b.PermutationTable[i]] = i;
    }
    for (sint16 i = 1; i <= GetParameter(); ++i)
        f.PermutationTable[i] = i;
    MeetSub(u, v, f.PermutationTable.data(), 1, GetParameter());

    return f;
};

bool ArtinBraidUnderlying::Compare(const ArtinBraidUnderlying &b) const {
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        if (PermutationTable[i] != b.PermutationTable[i]) {
            return false;
        }
    }
    return true;
};

ArtinBraidUnderlying ArtinBraidUnderlying::Inverse() const {
    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        f.PermutationTable[PermutationTable[i]] = i;
    }
    return f;
};

ArtinBraidUnderlying
ArtinBraidUnderlying::Product(const ArtinBraidUnderlying &b) const {
    ArtinBraidUnderlying f = ArtinBraidUnderlying(GetParameter());
    sint16 i;
    for (i = 1; i <= GetParameter(); i++) {
        f.PermutationTable[i] = b.PermutationTable[PermutationTable[i]];
    }
    return f;
};

ArtinBraidUnderlying
ArtinBraidUnderlying::LeftComplement(const ArtinBraidUnderlying &b) const {
    return b.Product(Inverse());
};

ArtinBraidUnderlying
ArtinBraidUnderlying::RightComplement(const ArtinBraidUnderlying &b) const {
    return Inverse().Product(b);
};

void ArtinBraidUnderlying::Randomize() {
    for (sint16 i = 1; i <= GetParameter(); ++i)
        PermutationTable[i] = i;
    for (sint16 i = 1; i < GetParameter(); ++i) {
        sint16 j = i + sint16(std::rand() / (RAND_MAX + 1.0) *
                              (GetParameter() - i + 1));
        sint16 z = PermutationTable[i];
        PermutationTable[i] = PermutationTable[j];
        PermutationTable[j] = z;
    }
};

std::vector<ArtinBraidUnderlying> ArtinBraidUnderlying::Atoms() const {
    sint16 i, n = GetParameter();
    ArtinBraidUnderlying atom(n);
    std::vector<ArtinBraidUnderlying> atoms;
    for (i = 1; i <= n - 1; i++) {
        atom.Identity();
        atom.PermutationTable[i] = i + 1;
        atom.PermutationTable[i + 1] = i;
        atoms.push_back(atom);
    }
    return atoms;
};

void ArtinBraidUnderlying::DeltaConjugate(sint16 k) {
    sint16 n = GetParameter();
    if (k % 2 != 0) {
        for (sint16 i = 1; i <= n / 2; i++) {
            sint16 u = PermutationTable[i];
            PermutationTable[i] = n - PermutationTable[n - i + 1] + 1;
            PermutationTable[n - i + 1] = n - u + 1;
        }
        if (n % 2 != 0) {
            PermutationTable[n / 2 + 1] = n - PermutationTable[n / 2 + 1] + 1;
        }
    }
}

} // namespace CGarside
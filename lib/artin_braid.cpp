#include "artin_braid.h"

namespace cgarside {

namespace artin {

Underlying::ParameterType Underlying::parameter_of_string(const std::string &str) {
    std::smatch match;

    if (std::regex_match(str, match,
                          std::regex{"[\\s\\t]*(" + number_regex + ")[\\s\\t]*"},
                          std::regex_constants::match_continuous)) {
        sint16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Number of strands is too big!\n" +
                                     match.str(1) + "could not be converted to a C++ integer.");
        }
        if (((2 <= i) && (i <= MaxBraidIndex))) {
            return i;
        } else if (2 > i) {
            throw InvalidStringError("Number of strands should be at least 2!");
        } else {
            throw InvalidStringError("Number of strands is too big!\n" +
                                     match.str(1) + " is strictly greater than " +
                                     std::to_string(MaxBraidIndex) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

void Underlying::MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s,
                         sint16 t) {
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

sint16 Underlying::GetParameter() const { return PresentationParameter; };

Underlying::Underlying(sint16 n)
    : PresentationParameter(n), PermutationTable(n + 1) {}

void Underlying::OfString(const std::string &str, size_t &pos) {
    sint16 n = GetParameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"(" + number_regex + ")"},
                          std::regex_constants::match_continuous)) {
        sint16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Invalid index for Artin generator!\n" +
                                     match.str(1) + " is not in [1, " +
                                     std::to_string(n) + "[.");
        }
        pos += match[0].length();
        if ((i >= 1) && (i < n)) {
            Identity();
            PermutationTable[i] = i + 1;
            PermutationTable[i + 1] = i;
        } else {
            throw InvalidStringError("Invalid index for Artin generator!\n" +
                                     match.str(1) + " is not in [1, " +
                                     std::to_string(n) + "[.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"D"},
                                 std::regex_constants::match_continuous)) {
        Delta();
    } else {
        throw InvalidStringError(
            std::string("Could not extract a factor from\n\"" + str.substr(pos) +
                        "\"!\nA factor should match regex ['1' - '9'] ['0' - '9']* | 'D'."));
    }
};

void Underlying::Print(IndentedOStream &os) const {
    sint16 i, j, k, n = GetParameter();

    Underlying c = Underlying(*this);

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

sint16 Underlying::LatticeHeight() const {
    sint16 n = GetParameter();
    return n * (n - 1) / 2;
}

void Underlying::Identity() {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = i;
    }
};

void Underlying::Delta() {
    sint16 i, n = GetParameter();
    for (i = 1; i <= n; i++) {
        PermutationTable[i] = n + 1 - i;
    }
};

Underlying Underlying::LeftMeet(const Underlying &b) const {
    thread_local sint16 s[MaxBraidIndex];

    Underlying f = Underlying(GetParameter());

    for (sint16 i = 1; i <= GetParameter(); ++i)
        s[i] = i;
    MeetSub(PermutationTable.data(), b.PermutationTable.data(), s, 1,
            GetParameter());
    for (sint16 i = 1; i <= GetParameter(); ++i)
        f.PermutationTable[s[i]] = i;

    return f;
};

Underlying Underlying::RightMeet(const Underlying &b) const {
    thread_local sint16 u[MaxBraidIndex], v[MaxBraidIndex];

    Underlying f = Underlying(GetParameter());

    for (sint16 i = 1; i <= GetParameter(); ++i) {
        u[PermutationTable[i]] = i;
        v[b.PermutationTable[i]] = i;
    }
    for (sint16 i = 1; i <= GetParameter(); ++i)
        f.PermutationTable[i] = i;
    MeetSub(u, v, f.PermutationTable.data(), 1, GetParameter());

    return f;
};

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

void Underlying::Randomize() {
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

std::vector<Underlying> Underlying::Atoms() const {
    sint16 i, n = GetParameter();
    Underlying atom(n);
    std::vector<Underlying> atoms;
    for (i = 1; i <= n - 1; i++) {
        atom.Identity();
        atom.PermutationTable[i] = i + 1;
        atom.PermutationTable[i + 1] = i;
        atoms.push_back(atom);
    }
    return atoms;
};

void Underlying::DeltaConjugate(sint16 k) {
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

size_t Underlying::Hash() const {
    size_t h = 0;
    for (sint16 i = 1; i <= GetParameter(); i++) {
        h = h * 31 + PermutationTable[i];
    }
    return h;
}

void Underlying::tableau(sint16 **&tab) const {
    sint16 i, j;
    sint16 n = GetParameter();
    for (i = 0; i < n; i++) {
        tab[i][i] = PermutationTable[i + 1];
    }
    for (j = 1; j <= n - 1; j++) {
        for (i = 0; i <= n - 1 - j; i++) {
            if (tab[i][i + j - 1] > tab[i + 1][i + j])
                tab[i][i + j] = tab[i][i + j - 1];
            else
                tab[i][i + j] = tab[i + 1][i + j];
        }
    }

    for (j = 1; j <= n - 1; j++) {
        for (i = j; i <= n - 1; i++) {
            if (tab[i - 1][i - j] < tab[i][i - j + 1])
                tab[i][i - j] = tab[i - 1][i - j];
            else
                tab[i][i - j] = tab[i][i - j + 1];
        }
    }
}

bool preserves_circles(const ArtinBraid &b) {
    sint16 j, k, t, d, n = b.GetParameter();
    sint16 *disj = new sint16[n + 1];

    Underlying delta_under(n);
    delta_under.Delta();

    sint16 cl = b.CanonicalLength();
    sint16 delta, itype = 0;
    if (b.Delta < 0)
        delta = -b.Delta;
    else
        delta = b.Delta;

    delta = delta % 2;

    sint16 ***tabarray = new sint16 **[cl + delta];
    ArtinBraid::ConstFactorItr it = b.FactorList.begin();

    for (j = 0; j < cl + delta; j++) {
        tabarray[j] = new sint16 *[n];
        for (k = 0; k < n; k++) {
            tabarray[j][k] = new sint16[n];
        }
        if (delta && j == 0)
            delta_under.tableau(tabarray[j]);
        else {
            (*it).GetUnderlying().tableau(tabarray[j]);
            it++;
        }
    }

    sint16 *bkmove = new sint16[n];
    sint16 bk;
    for (j = 2; j < n; j++) {
        for (k = 1; k <= n - j + 1; k++) {
            bk = k;
            for (t = 0; t < cl + delta; t++) {
                if (tabarray[t][bk - 1][j + bk - 2] -
                        tabarray[t][j + bk - 2][bk - 1] ==
                    j - 1)
                    bk = tabarray[t][j + bk - 2][bk - 1];
                else {
                    bk = 0;
                    break;
                }
            }
            if (bk == k) {
                itype = 1;
                j = n + 1;
                break;
            } else if (bk - k < j && k - bk < j)
                bk = 0;

            bkmove[k] = bk;
        }
        for (k = 1; k <= n - j + 1; k++) {
            for (d = 1; d <= n; d++)
                disj[d] = 1;

            bk = k;
            while (bk) {
                if (bkmove[bk] == k) {
                    itype = 1;
                    k = n - j;
                    j = n;
                    break;
                }
                for (d = bk - j + 1; d <= bk + j - 1; d++) {
                    if (d >= 1 && d <= n && d != k)
                        disj[d] = 0;
                }
                bk = bkmove[bk];
                if (disj[bk] == 0)
                    bk = 0;
            }
        }
    }

    if (itype)
        return true;
    else
        return false;
}

ThurstonType
thurston_type(const ArtinBraid &b,
              const ultra_summit::UltraSummitSet<ArtinBraid> &uss) {
    sint16 n = b.GetParameter();

    ArtinBraid pow = b;

    for (sint16 i = 0; i < n; i++) {
        if (pow.CanonicalLength() == 0)
            return ThurstonType::Periodic;
        pow.RightProduct(b);
    }

    for (typename ultra_summit::UltraSummitSet<ArtinBraid>::ConstIterator it =
             uss.begin();
         it != uss.end(); it++) {
        if (preserves_circles(*it)) {
            return ThurstonType::Reducible;
        }
    }

    return ThurstonType::PseudoAsonov;
}

ThurstonType thurston_type(const ArtinBraid &b) {
    return thurston_type(b, ultra_summit::USS(b));
}

} // namespace artin

template <>
IndentedOStream &IndentedOStream::operator<< <artin::ThurstonType>(
    const artin::ThurstonType &type) {
    switch (type) {
    case artin::ThurstonType::Periodic:
        os << "periodic";
        break;
    case artin::ThurstonType::Reducible:
        os << "reducible";
        break;
    case artin::ThurstonType::PseudoAsonov:
        os << "pseudo-Asonov";
        break;
    }
    return *this;
};

} // namespace cgarside
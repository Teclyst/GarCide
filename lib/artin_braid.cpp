#include "artin_braid.h"

namespace cgarside {

namespace artin {

Underlying::Parameter
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
        if (((2 <= i) && (i <= MAX_NUMBER_OF_STRANDS))) {
            return i;
        } else if (2 > i) {
            throw InvalidStringError("Number of strands should be at least 2!");
        } else {
            throw InvalidStringError("Number of strands is too big!\n" +
                                     match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MAX_NUMBER_OF_STRANDS) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

void Underlying::MeetSub(const sint16 *a, const sint16 *b, sint16 *r, sint16 s,
                         sint16 t) {
    thread_local sint16 u[MAX_NUMBER_OF_STRANDS], v[MAX_NUMBER_OF_STRANDS], w[MAX_NUMBER_OF_STRANDS];

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

sint16 Underlying::get_parameter() const { return number_of_strands; };

Underlying::Underlying(sint16 n)
    : number_of_strands(n), permutation_table(n + 1) {}

void Underlying::of_string(const std::string &str, size_t &pos) {
    Parameter n = get_parameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"(" + number_regex + ")"},
                          std::regex_constants::match_continuous)) {
        sint16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        pos += match[0].length();
        if ((i >= 1) && (i < n)) {
            identity();
            permutation_table[i] = i + 1;
            permutation_table[i + 1] = i;
        } else {
            throw InvalidStringError("Invalid index for Artin generator!\n" +
                                     match.str(1) + " is not in [1, " +
                                     std::to_string(n) + "[.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"D"},
                                 std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else {
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex Z | 'D',\nwhere Z matches "
            "integers."));
    }
};

void Underlying::print(IndentedOStream &os) const {
    sint16 i, j, k;
    Parameter n = get_parameter();

    Underlying c = Underlying(*this);

    bool is_first = true;

    for (i = 2; i <= n; i++) {
        for (j = i; j > 1 && c.permutation_table[j] < c.permutation_table[j - 1];
             j--) {
            os << (is_first ? "" : " ") << j - 1;
            is_first = false;
            k = c.permutation_table[j];
            c.permutation_table[j] = c.permutation_table[j - 1];
            c.permutation_table[j - 1] = k;
        }
    }
};

void Underlying::debug(IndentedOStream &os) const {
    os << "{   ";
    os.Indent(4);
    os << "number_of_strands:";
    os.Indent(4);
    os << EndLine() << get_parameter();
    os.Indent(-4);
    os << EndLine();
    os << "permutation_table:";
    os.Indent(4);
    os << EndLine();
    os << "[";
    for (sint16 i = 1; i < get_parameter(); i++) {
        os << permutation_table[i] << ", ";
    }
    os << permutation_table[get_parameter()];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

sint16 Underlying::lattice_height() const {
    Parameter n = get_parameter();
    return n * (n - 1) / 2;
}

void Underlying::identity() {
    Parameter n = get_parameter();
    for (sint16 i = 1; i <= n; i++) {
        permutation_table[i] = i;
    }
};

void Underlying::delta() {
    Parameter n = get_parameter();
    for (sint16 i = 1; i <= n; i++) {
        permutation_table[i] = n + 1 - i;
    }
};

Underlying Underlying::left_meet(const Underlying &b) const {
    thread_local sint16 s[MAX_NUMBER_OF_STRANDS];

    Underlying f = Underlying(get_parameter());

    for (sint16 i = 1; i <= get_parameter(); ++i)
        s[i] = i;
    MeetSub(permutation_table.data(), b.permutation_table.data(), s, 1,
            get_parameter());
    for (sint16 i = 1; i <= get_parameter(); ++i)
        f.permutation_table[s[i]] = i;

    return f;
};

Underlying Underlying::right_meet(const Underlying &b) const {
    thread_local sint16 u[MAX_NUMBER_OF_STRANDS], v[MAX_NUMBER_OF_STRANDS];

    Underlying f = Underlying(get_parameter());

    for (sint16 i = 1; i <= get_parameter(); ++i) {
        u[permutation_table[i]] = i;
        v[b.permutation_table[i]] = i;
    }
    for (sint16 i = 1; i <= get_parameter(); ++i)
        f.permutation_table[i] = i;
    MeetSub(u, v, f.permutation_table.data(), 1, get_parameter());

    return f;
};

bool Underlying::compare(const Underlying &b) const {
    sint16 i;
    for (i = 1; i <= get_parameter(); i++) {
        if (permutation_table[i] != b.permutation_table[i]) {
            return false;
        }
    }
    return true;
};

Underlying Underlying::inverse() const {
    Underlying f = Underlying(get_parameter());
    sint16 i;
    for (i = 1; i <= get_parameter(); i++) {
        f.permutation_table[permutation_table[i]] = i;
    }
    return f;
};

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    sint16 i;
    for (i = 1; i <= get_parameter(); i++) {
        f.permutation_table[i] = b.permutation_table[permutation_table[i]];
    }
    return f;
};

Underlying Underlying::left_complement(const Underlying &b) const {
    return b.product(inverse());
};

Underlying Underlying::right_complement(const Underlying &b) const {
    return inverse().product(b);
};

void Underlying::randomize() {
    for (sint16 i = 1; i <= get_parameter(); ++i)
        permutation_table[i] = i;
    for (sint16 i = 1; i < get_parameter(); ++i) {
        sint16 j = i + sint16(std::rand() / (RAND_MAX + 1.0) *
                              (get_parameter() - i + 1));
        sint16 z = permutation_table[i];
        permutation_table[i] = permutation_table[j];
        permutation_table[j] = z;
    }
};

std::vector<Underlying> Underlying::atoms() const {
    sint16 i;
    Parameter n = get_parameter();
    Underlying atom(n);
    std::vector<Underlying> atoms;
    for (i = 1; i <= n - 1; i++) {
        atom.identity();
        atom.permutation_table[i] = i + 1;
        atom.permutation_table[i + 1] = i;
        atoms.push_back(atom);
    }
    return atoms;
};

void Underlying::delta_conjugate_mut(sint16 k) {
    Parameter n = get_parameter();
    if (k % 2 != 0) {
        for (sint16 i = 1; i <= n / 2; i++) {
            sint16 u = permutation_table[i];
            permutation_table[i] = n - permutation_table[n - i + 1] + 1;
            permutation_table[n - i + 1] = n - u + 1;
        }
        if (n % 2 != 0) {
            permutation_table[n / 2 + 1] = n - permutation_table[n / 2 + 1] + 1;
        }
    }
}

size_t Underlying::hash() const {
    size_t h = 0;
    for (sint16 i = 1; i <= get_parameter(); i++) {
        h = h * 31 + permutation_table[i];
    }
    return h;
}

void Underlying::tableau(sint16 **&tab) const {
    sint16 i, j;
    Braid::Parameter n = get_parameter();
    for (i = 0; i < n; i++) {
        tab[i][i] = permutation_table[i + 1];
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

bool preserves_circles(const Braid &b) {
    sint16 j, k, t, d;
    Braid::Parameter n = b.get_parameter();
    sint16 *disj = new sint16[n + 1];

    Underlying delta_under(n);
    delta_under.delta();

    sint16 cl = int(b.canonical_length());
    sint16 delta, itype = 0;
    if (b.inf() < 0)
        delta = -b.inf();
    else
        delta = b.inf();

    delta = delta % 2;

    sint16 ***tabarray = new sint16 **[cl + delta];
    Braid::ConstFactorItr it = b.cbegin();

    for (j = 0; j < cl + delta; j++) {
        tabarray[j] = new sint16 *[n];
        for (k = 0; k < n; k++) {
            tabarray[j][k] = new sint16[n];
        }
        if (delta && j == 0)
            delta_under.tableau(tabarray[j]);
        else {
            (*it).get_underlying().tableau(tabarray[j]);
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
thurston_type(const Braid &b,
              const ultra_summit::UltraSummitSet<Braid> &uss) {
    Braid::Parameter n = b.get_parameter();

    Braid pow = b;

    for (sint16 i = 0; i < n; i++) {
        if (pow.canonical_length() == 0)
            return ThurstonType::Periodic;
        pow.right_multiply(b);
    }

    for (typename ultra_summit::UltraSummitSet<Braid>::ConstIterator it =
             uss.begin();
         it != uss.end(); it++) {
        if (preserves_circles(*it)) {
            return ThurstonType::Reducible;
        }
    }

    return ThurstonType::PseudoAsonov;
}

ThurstonType thurston_type(const Braid &b) {
    return thurston_type(b, ultra_summit::ultra_summit_set(b));
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
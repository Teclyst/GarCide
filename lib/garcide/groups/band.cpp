/**
 * @file band.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for standard braid groups (dual Garside
 * structure).
 * @version 0.1
 * @date 2024-07-28
 *
 * @copyright Copyright (C) 2024. Distributed under the GNU General Public
 * License, version 3.
 *
 */

/*
 * GarCide Copyright (C) 2024 Matteo Wei.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License in LICENSE for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "garcide/groups/band.hpp"

namespace garcide::band {

Underlying::Parameter Underlying::parameter_of_string(const std::string &str) {
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
            throw InvalidStringError(
                "Number of strands is too big!\n" + match.str(1) +
                " is strictly greater than " +
                std::to_string(MAX_NUMBER_OF_STRANDS) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
}

Underlying::Underlying(sint16 n) : permutation_table(n + 1) {}

void Underlying::print(IndentedOStream &os) const {
    // Recall that a band braid is represented by decreasing cycles.
    sint16 i, j, n = get_parameter();
    std::vector<sint16> curr_cycle;
    std::vector<bool> seen(n + 1, false);
    bool is_first = true;
    for (i = 1; i <= n; ++i) {
        if (!seen[i]) {
            j = i;
            curr_cycle.clear();
            while (j < permutation_table[j]) {
                curr_cycle.push_back(j);
                seen[j] = true;
                j = permutation_table[j];
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

void Underlying::of_string(const std::string &str, size_t &pos) {
    sint16 n = get_parameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
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
            identity();
            permutation_table[i] = j;
            permutation_table[j] = i;
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

void Underlying::assign_partition(sint16 *x) const {
    for (sint16 i = 1; i <= get_parameter(); ++i)
        x[i] = 0;
    for (sint16 i = 1; i <= get_parameter(); ++i) {
        if (x[i] == 0)
            x[i] = i;
        if (permutation_table[i] > i)
            x[permutation_table[i]] = x[i];
    }
}

void Underlying::of_partition(const sint16 *x) {
    thread_local sint16 z[MAX_NUMBER_OF_STRANDS];

    for (sint16 i = 1; i <= get_parameter(); ++i)
        z[i] = 0;
    for (sint16 i = get_parameter(); i >= 1; --i) {
        permutation_table[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
        z[x[i]] = i;
    }
}

Underlying Underlying::left_meet(const Underlying &b) const {
    thread_local sint16 x[MAX_NUMBER_OF_STRANDS], y[MAX_NUMBER_OF_STRANDS],
        z[MAX_NUMBER_OF_STRANDS];

    assign_partition(x);
    b.assign_partition(y);

    thread_local sint16 P[MAX_NUMBER_OF_STRANDS][MAX_NUMBER_OF_STRANDS];

    for (sint16 i = get_parameter(); i >= 1; i--) {
        P[x[i]][y[i]] = i;
    }

    for (sint16 i = 1; i <= get_parameter(); i++) {
        z[i] = P[x[i]][y[i]];
    }

    Underlying c = Underlying(*this);

    c.of_partition(z);

    return c;
}

void Underlying::identity() {
    sint16 i, n = get_parameter();
    for (i = 1; i <= n; i++) {
        permutation_table[i] = i;
    }
}

void Underlying::delta() {
    sint16 i, n = get_parameter();
    for (i = 1; i < n; i++) {
        permutation_table[i] = i + 1;
    }
    permutation_table[n] = 1;
}

Underlying Underlying::inverse() const {
    Underlying f = Underlying(get_parameter());
    sint16 i;
    for (i = 1; i <= get_parameter(); i++) {
        f.permutation_table[permutation_table[i]] = i;
    }
    return f;
}

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    sint16 i;
    for (i = 1; i <= get_parameter(); i++) {
        f.permutation_table[i] = b.permutation_table[permutation_table[i]];
    }
    return f;
}

Underlying Underlying::left_complement(const Underlying &b) const {
    return b.product(inverse());
}

Underlying Underlying::right_complement(const Underlying &b) const {
    return inverse().product(b);
}

void Underlying::delta_conjugate_mut(sint16 k) {
    Underlying under = *this;
    sint16 i, n = get_parameter();

    for (i = 1; i <= n; i++) {
        under.permutation_table[i] =
            Rem(permutation_table[Rem(i - k - 1, n) + 1] + k - 1, n) + 1;
    }
    *this = under;
}

void Underlying::randomize() {
#ifdef USE_CLN
    static sint8 s[2 * MAX_NUMBER_OF_STRANDS + 1];
    cln::cl_I k = cln::random_I(cln::default_random_state,
                                get_catalan_number(get_parameter())) +
                  1;
    ballot_sequence(get_parameter(), k, s);
    of_ballot_sequence(s);
#else
    throw NonRandomizable();
#endif
}

std::vector<Underlying> Underlying::atoms() const {
    sint16 n = get_parameter();
    Underlying atom(n);
    atom.identity();
    std::vector<Underlying> atoms;
    for (sint16 i = 1; i <= n; i++) {
        for (sint16 j = 1; j < i; j++) {
            atom.permutation_table[i] = j;
            atom.permutation_table[j] = i;
            atoms.push_back(atom);
            atom.permutation_table[i] = i;
            atom.permutation_table[j] = j;
        }
    }
    return atoms;
}

size_t Underlying::hash() const {
    std::size_t h = 0;
    for (sint16 i = 1; i <= get_parameter(); i++) {
        h = h * 31 + permutation_table[i];
    }
    return h;
}

void Underlying::of_ballot_sequence(const sint8 *s) {
    static sint16 stack[MAX_NUMBER_OF_STRANDS];
    sint16 sp = 0;

    for (sint16 i = 1; i <= 2 * get_parameter(); ++i) {
        if (s[i] == 1) {
            stack[sp++] = i;
        } else {
            sint16 j = stack[--sp];
            if ((i / 2) * 2 != i)
                permutation_table[j / 2] = (i + 1) / 2;
            else
                permutation_table[i / 2] = (j + 1) / 2;
        }
    }
}

#ifdef USE_CLN

void ballot_sequence(sint16 n, cln::cl_I k, sint8 *s) {
    sint16 i;
    cln::cl_I r;

    if (k <= (r = get_catalan_number(n - 1) * get_catalan_number(0)))
        i = 1;
    else if (k > (r = get_catalan_number(n) - r)) {
        i = n;
        k = k - r;
    } else {
        for (i = 1; i <= n; ++i) {
            if (k <=
                (r = get_catalan_number(i - 1) * get_catalan_number(n - i)))
                break;
            else
                k = k - r;
        }
    }

    cln::cl_I_div_t d = cln::floor2(k - 1, get_catalan_number(n - i));

    s[1] = 1;
    s[2 * i] = -1;
    if (i > 1)
        ballot_sequence(i - 1, d.quotient + 1, s + 1);
    if (i < n)
        ballot_sequence(n - i, d.remainder + 1, s + 2 * i);
}

class _CatalanNumber {
    friend const cln::cl_I &get_catalan_number(sint16 n);

  public:
    _CatalanNumber();

  private:
    cln::cl_I table[Underlying::MAX_NUMBER_OF_STRANDS + 1];
    cln::cl_I &c(sint16 n) { return table[n]; }
}

static _CatalanNumber CatalanNumber;

_CatalanNumber::_CatalanNumber() {
    c(0) = 1;
    for (sint16 n = 1; n < Underlying::MAX_NUMBER_OF_STRANDS; ++n) {
        c(n) = 0;
        for (sint16 k = 0; k < n; ++k) {
            c(n) = c(n) + c(k) * c(n - k - 1);
        }
    }
}

const cln::cl_I &get_catalan_number(sint16 n) { return CatalanNumber.c(n); }

#endif

} // namespace garcide::band

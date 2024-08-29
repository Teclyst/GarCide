/**
 * @file octahedral.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for \f$\mathbf B\f$-series Artin groups (dual
 * Garside structure).
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

#include "garcide/groups/octahedral.hpp"

namespace garcide::octahedral {

Underlying::Parameter Underlying::parameter_of_string(const std::string &str) {
    std::smatch match;

    if (std::regex_match(str, match,
                         std::regex{"[\\s\\t]*(" + number_regex + ")[\\s\\t]*"},
                         std::regex_constants::match_continuous)) {
        i16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Parameter is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        if (((1 <= i) && (i <= MAX_PARAMETER))) {
            return i;
        } else if (1 > i) {
            throw InvalidStringError("Parameter should be at least 1!");
        } else {
            throw InvalidStringError("Parameter strands is too big!\n" +
                                     match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MAX_PARAMETER) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
}

Underlying::Underlying(i16 n) : permutation_table(2 * n + 1) {}

void Underlying::print(IndentedOStream &os) const {
    // Recall that a band braid is represented by decreasing cycles.
    i16 i, j, n = get_parameter();
    std::vector<i16> curr_cycle;
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
                j = permutation_table[j];
            }
            if (!is_first && int(curr_cycle.size()) > 1) {
                os << " ";
            } else if (int(curr_cycle.size()) > 1) {
                is_first = false;
            }
            for (i16 l = int(curr_cycle.size()) - 1; l >= 1; --l) {
                os << "s(" << curr_cycle[l] << ", " << curr_cycle[l - 1] << ")"
                   << ((l == 1) ? "" : " ");
            }
            // Long cycle.
            if (j > n) {
                os << (is_first ? "l" : " l") << j - n;
            }
        }
    }
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
    for (i16 i = 1; i < 2 * get_parameter(); i++) {
        os << permutation_table[i] << ", ";
    }
    os << permutation_table[2 * get_parameter()];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

void Underlying::of_string(const std::string &str, size_t &pos) {
    i16 n = get_parameter();
    std::smatch match;
    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else if (std::regex_search(
                   str.begin() + pos, str.end(), match,
                   std::regex{"(?:s[\\s\\t]*_?[\\s\\t]*)?\\([\\s\\t]*(" +
                              number_regex + ")[\\s\\t]*,?[\\s\\t]*(" +
                              number_regex + ")[\\s\\t]*\\)"},
                   std::regex_constants::match_continuous)) {
        i16 i, j;
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
            identity();
            permutation_table[i] = j;
            permutation_table[j] = i;
            permutation_table[Rem(i + n - 1, 2 * n) + 1] =
                Rem(j + n - 1, 2 * n) + 1;
            permutation_table[Rem(j + n - 1, 2 * n) + 1] =
                Rem(i + n - 1, 2 * n) + 1;
        } else {
            throw InvalidStringError(
                "Indexes for short generators should not be equal mod " +
                std::to_string(n) + "!\n(" + match.str(1) + ", " +
                match.str(2) + ") is not a valid factor.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"(?:l[\\s\\t]*_?[\\s\\t]*)?(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        i16 i;
        try {
            i = std::stoi(match[1]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        i = Rem(i - 1, 2 * n) + 1;
        pos += match[0].length();
        identity();
        permutation_table[i] = Rem(i + n - 1, 2 * n) + 1;
        permutation_table[Rem(i + n - 1, 2 * n) + 1] = i;
    } else {
        throw InvalidStringError(
            "Could not extract a factor from \"" + str.substr(pos) +
            "\"!\nA factor should match regex\n('s' '_'?)? '(' Z ','? Z ')' | "
            "('l' '_'?)? Z "
            "| "
            "'D',\nwhere Z matches integers, and ignoring whitespaces.");
    }
}

void Underlying::assign_partition(i16 *x) const {
    for (i16 i = 1; i <= 2 * get_parameter(); ++i)
        x[i] = 0;
    for (i16 i = 1; i <= 2 * get_parameter(); ++i) {
        if (x[i] == 0)
            x[i] = i;
        if (permutation_table[i] > i)
            x[permutation_table[i]] = x[i];
    }
}

void Underlying::of_partition(const i16 *x) {
    thread_local i16 z[2 * MAX_PARAMETER + 1];

    for (i16 i = 1; i <= 2 * get_parameter(); ++i)
        z[i] = 0;
    for (i16 i = 2 * get_parameter(); i >= 1; --i) {
        permutation_table[i] = (z[x[i]] == 0) ? x[i] : z[x[i]];
        z[x[i]] = i;
    }
}

Underlying Underlying::left_meet(const Underlying &b) const {
    thread_local i16 x[2 * MAX_PARAMETER + 1], y[2 * MAX_PARAMETER + 1],
        z[2 * MAX_PARAMETER + 1];

    assign_partition(x);
    b.assign_partition(y);

    thread_local i16 P[2 * MAX_PARAMETER + 1][2 * MAX_PARAMETER + 1];

    for (i16 i = 2 * get_parameter(); i >= 1; i--) {
        P[x[i]][y[i]] = i;
    }

    for (i16 i = 1; i <= 2 * get_parameter(); i++) {
        z[i] = P[x[i]][y[i]];
    }

    Underlying c = Underlying(*this);

    c.of_partition(z);

    return c;
}

void Underlying::identity() {
    for (i16 i = 1; i <= 2 * get_parameter(); i++) {
        permutation_table[i] = i;
    }
}

void Underlying::delta() {
    i16 i, n = get_parameter();
    for (i = 1; i < 2 * n; i++) {
        permutation_table[i] = i + 1;
    }
    permutation_table[2 * n] = 1;
}

Underlying Underlying::inverse() const {
    Underlying f = Underlying(get_parameter());
    i16 i;
    for (i = 1; i <= 2 * get_parameter(); i++) {
        f.permutation_table[permutation_table[i]] = i;
    }
    return f;
}

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    i16 i;
    for (i = 1; i <= 2 * get_parameter(); i++) {
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

void Underlying::delta_conjugate_mut(i16 k) {
    Underlying under = *this;
    i16 i, n = get_parameter();

    for (i = 1; i <= 2 * n; i++) {
        under.permutation_table[i] =
            Rem(permutation_table[Rem(i - k - 1, 2 * n) + 1] + k - 1, 2 * n) +
            1;
    }
    *this = under;
}

std::vector<Underlying> Underlying::atoms() const {
    i16 n = get_parameter();
    std::vector<Underlying> atoms;
    for (i16 i = 1; i <= n; i++) {
        for (i16 j = i + 1; j <= n; j++) {
            Underlying atom = Underlying(n);
            atom.identity();
            atom.permutation_table[i] = j;
            atom.permutation_table[j] = i;
            atom.permutation_table[i + n] = j + n;
            atom.permutation_table[j + n] = i + n;
            atoms.push_back(atom);
        }
        for (i16 j = n + i + 1; j <= 2 * n; j++) {
            Underlying atom = Underlying(n);
            atom.identity();
            atom.permutation_table[i] = j;
            atom.permutation_table[j] = i;
            atom.permutation_table[i + n] = j - n;
            atom.permutation_table[j - n] = i + n;
            atoms.push_back(atom);
        }
        Underlying atom = Underlying(n);
        atom.identity();
        atom.permutation_table[i] = n + i;
        atom.permutation_table[n + i] = i;
        atoms.push_back(atom);
    }
    return atoms;
}

std::size_t Underlying::hash() const {
    std::size_t h = 0;
    for (i16 i = 1; i <= 2 * get_parameter(); i++) {
        h = h * 31 + permutation_table[i];
    }
    return h;
}

} // namespace garcide::octahedral
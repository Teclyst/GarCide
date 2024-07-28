/**
 * @file dual_complex.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for B(e, e, n) groups (dual Garside structure).
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

#include "dual_complex.h"

namespace garcide {

namespace dual_complex {

void EENParameter::print(IndentedOStream &os) const {
    os << "(e: " << e << ", n: " << n << ")";
}

EENParameter Underlying::get_parameter() const { return een_index; }

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
        if ((2 <= e) && (2 <= n) && (n <= MAX_N) && (e * n <= MAX_E * MAX_N)) {
            return EENParameter(e, n);
        } else if (2 > e) {
            throw InvalidStringError("e should be at least 2!");
            } else if (2 > n) {
            throw InvalidStringError("n should be at least 2!");
        } else if (n > MAX_N) {
            throw InvalidStringError("n is too big!\n" +
                                     match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MAX_N) + ".");
        } else {
            throw InvalidStringError("e * n is too big!\n" +
                                     match.str(1) + " * " + match.str(2) + " = " + std::to_string(e * n) +
                                     " is strictly greater than " +
                                     std::to_string(MAX_N * MAX_E) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
}

sint16 Underlying::lattice_height() const { return get_parameter().n + 1; }

Underlying::Underlying(Parameter p)
    : een_index(p), permutation_table(p.n + 1),
      coefficient_table(p.n + 1) {}

void Underlying::debug(IndentedOStream &os) const {
    os << "{   ";
    os.Indent(4);
    os << "een_index:";
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
        os << coefficient_table[i] << ", ";
    }
    os << coefficient_table[get_parameter().n];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

void Underlying::print(IndentedOStream &os) const {
    sint16 n = get_parameter().n, e = get_parameter().e;
    std::vector<bool> seen(n + 1, false);
    std::vector<sint16> curr_cycle;
    sint16 other_smallest = 0, cycle_type, curr, c = 0;
    seen[0] = true;
    bool is_first = true;

    curr = permutation_table[0];
    // Short assymetric case.
    if (curr != 0) {
        while (curr != 0) {
            curr_cycle.push_back(curr);
            seen[curr] = true;
            other_smallest =
                ((curr < curr_cycle[other_smallest]) && (curr != 0))
                    ? c
                    : other_smallest;
            curr = permutation_table[curr];
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
                       ? curr_cycle[i] + Rem(coefficient_table[0] + 1, e) * n
                       : curr_cycle[i] + coefficient_table[0] * n)
               << ", "
               << ((i >= other_smallest + 1)
                       ? curr_cycle[i - 1] + Rem(coefficient_table[0] + 1, e) * n
                       : curr_cycle[i - 1] + coefficient_table[0] * n)
               << ") ";
        }
        os << curr_cycle[0] + coefficient_table[0] * n;
    }

    for (sint16 i = 1; i <= n; ++i) {
        if (!seen[i]) {
            seen[i] = true;
            c = 0;
            other_smallest = 0;
            cycle_type = coefficient_table[i];
            if (coefficient_table[i] == e - 1) {
                other_smallest = c + 1;
            }
            curr = permutation_table[i];
            curr_cycle.clear();
            curr_cycle.push_back(i);
            while (curr != i) {
                curr_cycle.push_back(curr);
                seen[curr] = true;
                c++;
                if (coefficient_table[curr] == e - 1) {
                    other_smallest = c + 1;
                }
                cycle_type += coefficient_table[curr];
                curr = permutation_table[curr];
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

void Underlying::of_string(const std::string &str, size_t &pos) {
    sint16 n = get_parameter().n, e = get_parameter().e;
    std::smatch match;
    if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"D"},
                                 std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
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
            identity();
            permutation_table[i] = (j > n) ? j - n : j;
            permutation_table[(j > n) ? j - n : j] = i;
            coefficient_table[i] = (j > n) ? 1 : 0;
            coefficient_table[(j > n) ? j - n : j] = (j > n) ? e - 1 : 0;
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
        identity();
        permutation_table[0] = r;
        permutation_table[r] = 0;
        coefficient_table[0] = q;
        coefficient_table[r] = (q == 0) ? 0 : e - q;
    } else {
        throw InvalidStringError(
            "Could not extract a factor from \"" + str.substr(pos) +
            "\"!\nA factor should match regex '(' Z ','? Z ')' | Z | "
            "'D',\nwhere Z matches integers, and ignoring whitespaces.");
    }
}

void Underlying::assign_partition(sint16 *x) const {
    sint16 n = get_parameter().n, e = get_parameter().e;
    std::vector<sint16> curr_cycle;
    sint16 other_smallest, cycle_type, curr, c;
    for (sint16 i = 1; i <= e * n; ++i) {
        x[i] = -1;
    }
    x[0] = 0;

    curr = permutation_table[0];
    if (curr != 0) {
        other_smallest = 0;
        c = 0;
        while (curr != 0) {
            curr_cycle.push_back(curr);
            other_smallest =
                ((curr < curr_cycle[other_smallest]) && (curr != 0))
                    ? c
                    : other_smallest;
            curr = permutation_table[curr];
            c++;
        }
        if (other_smallest != 0) {
            for (sint16 k = 0; k < other_smallest; k++) {
                for (sint16 l = 0; l < e - 1; l++) {
                    x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
                }
                x[curr_cycle[k] + (e - 1) * n] = curr_cycle[other_smallest];
                x[curr_cycle[k] + Rem(coefficient_table[0], e) * n] = 0;
            }
            for (sint16 k = other_smallest; k < int(curr_cycle.size()); k++) {
                for (sint16 l = 0; l < e - 1; l++) {
                    x[curr_cycle[k] + (l + 1) * n] = curr_cycle[0] + l * n;
                }
                x[curr_cycle[k]] = curr_cycle[other_smallest];
                x[curr_cycle[k] + Rem(coefficient_table[0] + 1, e) * n] = 0;
            }
        } else {
            for (sint16 k = 0; k < int(curr_cycle.size()); k++) {
                for (sint16 l = 0; l < e; l++) {
                    x[curr_cycle[k] + l * n] = curr_cycle[0] + l * n;
                }
                x[curr_cycle[k] + Rem(coefficient_table[0], e) * n] = 0;
            }
        }
    }

    for (sint16 i = 1; i <= n; ++i) {
        if (x[i] < 0) {
            c = 0;
            other_smallest = 0;
            cycle_type = coefficient_table[i];
            if (coefficient_table[i] == e - 1) {
                other_smallest = 1;
            }
            curr = permutation_table[i];
            curr_cycle.clear();
            curr_cycle.push_back(i);
            while (curr != i) {
                curr_cycle.push_back(curr);
                c++;
                if (coefficient_table[curr] == e - 1) {
                    other_smallest = c + 1;
                }
                cycle_type += coefficient_table[curr];
                curr = permutation_table[curr];
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
void Underlying::of_partition(const sint16 *x) {
    thread_local sint16 z[MAX_N + 1];
    sint16 min_cycle_0 = 0, max_cycle_0 = 0;
    sint16 n = get_parameter().n, e = get_parameter().e, r;

    for (sint16 i = 0; i <= n; ++i) {
        z[i] = -1;
        permutation_table[i] = -1;
        coefficient_table[i] = -1;
    }
    // First find short symmetrical cycles.
    for (sint16 i = 2 * n - 1; i >= n + 1; --i) {
        if ((x[i] <= n) && (x[i] >= 1)) {
            r = i - n;
            if (z[x[i]] == -1) {
                permutation_table[r] = x[i];
                coefficient_table[r] = e - 1;
                z[x[i]] = r;
            } else {
                permutation_table[r] = z[x[i]];
                coefficient_table[r] = 0;
                z[x[i]] = r;
            }
        }
    }
    for (sint16 i = n; i >= 1; --i) {
        if ((x[i] <= n) && (x[i] >= 1) && (x[i + n] > n)) {
            if ((z[x[i]] == -1)) {
                permutation_table[i] = x[i];
                coefficient_table[i] = 0;
            } else {
                coefficient_table[i] = (z[x[i]] < i) ? 1 : 0;
                permutation_table[i] = z[x[i]];
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
            coefficient_table[0] = e - 1;
            permutation_table[0] = 0;
            for (sint16 i = n; i >= 1; i--) {
                if (x[i] == 0) {
                    if (z[x[i]] == -1) {
                        permutation_table[i] = min_cycle_0;
                        coefficient_table[i] = 1;
                    } else {
                        permutation_table[i] = z[x[i]];
                        coefficient_table[i] = 0;
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
            coefficient_table[0] = q_min;
            permutation_table[0] = r_min;

            for (sint16 i = n - 1; i >= n - r_min + 1; --i) {
                sint16 i_en = Rem(i + min_cycle_0 - 1, e * n) + 1;
                r = Rem(i + min_cycle_0 - 1, n) + 1;
                if (x[i_en] == 0) {
                    permutation_table[r] = z[0];
                    coefficient_table[r] = (z[0] == 0) ? Rem(e - q_max, e) : 0;
                    z[0] = r;
                }
            }
            for (sint16 i = n - r_min; i >= 0; --i) {
                sint16 i_en = Rem(i + min_cycle_0 - 1, e * n) + 1;
                r = Rem(i + min_cycle_0 - 1, n) + 1;
                if ((x[i_en] == 0)) {
                    if (z[0] == 0) {
                        permutation_table[r] = z[0];
                        coefficient_table[r] = Rem(e - q_min, e);
                    } else {
                        coefficient_table[r] = (z[0] < r) ? 1 : 0;
                        permutation_table[r] = z[0];
                    }
                    z[0] = r;
                }
            }
        }
    }

    // Short symmetric case.
    else {
        coefficient_table[0] = 0;
        permutation_table[0] = 0;
    }
}

Underlying Underlying::left_meet(const Underlying &b) const {
    thread_local sint16 x[MAX_E * MAX_N + 1], y[MAX_E * MAX_N + 1],
        z[MAX_E * MAX_N + 1];

    assign_partition(x);
    b.assign_partition(y);

    thread_local sint16 P[MAX_E * MAX_N + 1][MAX_E * MAX_N + 1];

    for (sint16 i = get_parameter().e * get_parameter().n; i >= 0; i--) {
        P[x[i]][y[i]] = i;
    }

    for (sint16 i = 0; i <= get_parameter().e * get_parameter().n; i++) {
        z[i] = P[x[i]][y[i]];
    }

    Underlying c = Underlying(*this);

    c.of_partition(z);

    return c;
}

void Underlying::identity() {
    for (sint16 i = 0; i <= get_parameter().n; i++) {
        permutation_table[i] = i;
        coefficient_table[i] = 0;
    }
}

void Underlying::delta() {
    sint16 i, n = get_parameter().n;
    for (i = 1; i <= n; i++) {
        permutation_table[i] = i + 1;
        coefficient_table[i] = 0;
    }
    permutation_table[0] = 0;
    coefficient_table[0] = get_parameter().e - 1;
    permutation_table[n] = 1;
    coefficient_table[n] = 1;
}

bool Underlying::compare(const Underlying &b) const {
    sint16 i;
    for (i = 0; i <= get_parameter().n; i++) {
        if ((permutation_table[i] != b.permutation_table[i]) ||
            (coefficient_table[i] != b.coefficient_table[i])) {
            return false;
        }
    }
    return true;
};

Underlying Underlying::inverse() const {
    Underlying f = Underlying(get_parameter());
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 0; i <= n; i++) {
        if ((permutation_table[i] > n) || (permutation_table[i] < 0)) {
            while (true) {
            };
        }
        f.permutation_table[permutation_table[i]] = i;
        f.coefficient_table[permutation_table[i]] =
            coefficient_table[i] == 0 ? 0 : e - coefficient_table[i];
    }
    return f;
};

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 0; i <= n; i++) {
        f.permutation_table[i] = b.permutation_table[permutation_table[i]];
        f.coefficient_table[i] = Rem(
            b.coefficient_table[permutation_table[i]] + coefficient_table[i], e);
    }
    return f;
};

Underlying Underlying::left_complement(const Underlying &b) const {
    return b.product(inverse());
};

Underlying Underlying::right_complement(const Underlying &b) const {
    return inverse().product(b);
};

void Underlying::delta_conjugate_mut(sint16 k) {
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    sint16 q = Quot(k, n), r = Rem(k, n), q_e = Rem(q, e);

    Underlying delta_k = Underlying(get_parameter());

    delta_k.permutation_table[0] = 0;
    delta_k.coefficient_table[0] = Rem(-k, e);
    for (i = 1; i <= n - r; i++) {
        delta_k.permutation_table[i] = Rem(i + k - 1, n) + 1;
        delta_k.coefficient_table[i] = q_e;
    }
    q_e += 1;
    q_e = ((q_e == e) ? 0 : q_e);
    for (i = n - r + 1; i <= n; i++) {
        delta_k.permutation_table[i] = Rem(i + k - 1, n) + 1;
        delta_k.coefficient_table[i] = q_e;
    }

    *this = delta_k.inverse().product((*this).product(delta_k));
}

void Underlying::randomize() { throw NonRandomizable(); }

std::vector<Underlying> Underlying::atoms() const {
    Parameter p = get_parameter();
    std::vector<Underlying> atoms;
    Underlying atom = Underlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 1; i <= n; i++) {
        for (sint16 j = i + 1; j <= n; j++) {
            atom.identity();
            atom.permutation_table[i] = j;
            atom.permutation_table[j] = i;
            atoms.push_back(atom);
        }
        for (sint16 j = 1; j < i; j++) {
            atom.identity();
            atom.permutation_table[i] = j;
            atom.permutation_table[j] = i;
            atom.coefficient_table[i] = 1;
            atom.coefficient_table[j] = e - 1;
            atoms.push_back(atom);
        }
    }
    for (sint16 i = 1; i <= n; i++) {
        atom.identity();
        atom.permutation_table[0] = i;
        atom.permutation_table[i] = 0;
        atom.coefficient_table[0] = 0;
        atom.coefficient_table[i] = 0;
        atoms.push_back(atom);
    }
    for (sint16 k = 1; k < e; k++) {
        for (sint16 i = 1; i <= n; i++) {
            atom.identity();
            atom.permutation_table[0] = i;
            atom.permutation_table[i] = 0;
            atom.coefficient_table[0] = k;
            atom.coefficient_table[i] = e - k;
            atoms.push_back(atom);
        }
    }
    return atoms;
}

std::size_t Underlying::hash() const {
    std::size_t h = 0;
    for (sint16 i = 1; i <= get_parameter().n; i++) {
        h = h * 31 + permutation_table[i];
    }
    for (sint16 i = 1; i <= get_parameter().n; i++) {
        h = h * 31 + coefficient_table[i];
    }
    return h;
}

} // namespace dual_complex

template <>
IndentedOStream &
    IndentedOStream::operator<< <dual_complex::EENParameter>(const dual_complex::EENParameter &p) {
    p.print(*this);
    return *this;
};

} // namespace cgarside
#include "standard_complex_reflection.h"
#include <iostream>
#include <string>

namespace cgarside {

namespace standard_complex {

void EENParameter::print(IndentedOStream &os) const {
    os << "(e: " << e << ", n: " << n << ")";
}

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
    os << "coefficient_table:";
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
        if ((2 <= e) && (2 <= n) && (n <= MAX_N)) {
            return EENParameter(e, n);
        } else if (2 > e) {
            throw InvalidStringError("e should be at least 2!");
        } else if (2 > n) {
            throw InvalidStringError("n should be at least 2!");
        } else {
            throw InvalidStringError("n is too big!\n" + match.str(1) +
                                     " is strictly greater than " +
                                     std::to_string(MAX_N) + ".");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

sint16 Underlying::lattice_height() const {
    sint16 n = get_parameter().n;
    return n * (n + 1);
}

Underlying::Underlying(Parameter p)
    : een_index(p), permutation_table(p.n), coefficient_table(p.n) {}

void Underlying::print(IndentedOStream &os) const {
    thread_local sint16 dir_perm[MAX_N];
    Underlying copy = *this;
    copy.direct(dir_perm);
    sint16 n = get_parameter().n, e = get_parameter().e;
    bool is_first = true;
    for (sint16 i = 2; i <= n; i++) {
        for (sint16 j = i; j > 2; j--) {
            if (copy.is_s_left_divisor(dir_perm, j)) {
                if (is_first) {
                    is_first = false;
                } else {
                    os << " ";
                }
                copy.s_left_multiply(dir_perm, j);
                os << "s" << j;
            } else {
                goto end_subword;
            }
        }

        if (copy.is_t_left_divisor(dir_perm, 1)) {
            if (is_first) {
                is_first = false;
            } else {
                os << " ";
            }
            copy.t_left_multiply(dir_perm, 1);
            os << "t" << 1;

            if (copy.is_t_left_divisor(dir_perm, 0)) {
                copy.t_left_multiply(dir_perm, 0);
                os << " t" << 0;

                for (sint16 j = 3; j < i + 1; j++) {
                    if (copy.is_s_left_divisor(dir_perm, j)) {
                        copy.s_left_multiply(dir_perm, j);
                        os << " s" << j;
                    } else {
                        goto end_subword;
                    }
                }
            }
            goto end_subword;
        }

        for (sint16 k = 0; k < e; k++) {
            if (copy.is_t_left_divisor(dir_perm, k)) {
                if (is_first) {
                    is_first = false;
                } else {
                    os << " ";
                }
                copy.t_left_multiply(dir_perm, k);
                os << "t" << k;
                goto end_subword;
            }
        }
    end_subword:;
    }
}

void Underlying::of_string(const std::string &str, size_t &pos) {
    sint16 n = get_parameter().n, e = get_parameter().e;
    std::smatch match;
    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"([st])[\\s\\t]*_?[\\s\\t]*(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        sint16 i = std::stoi(match[2]);
        pos += match[0].length();
        if (match[1] == "s") {
            if ((3 <= i) && (i <= n)) {
                // s_i.
                identity();
                permutation_table[i - 2] = i - 1;
                permutation_table[i - 1] = i - 2;
            } else {
                throw InvalidStringError(
                    "Invalid index for s type generator!\n" +
                    std::to_string(i) + " is not in [3, " + std::to_string(n) +
                    "].");
            }
        } else {
            // t_i.
            i = Rem(i, e);
            identity();
            permutation_table[0] = 1;
            permutation_table[1] = 0;
            coefficient_table[1] = i;
            coefficient_table[0] = e - ((i == 0) ? e : i);
        }
    } else {
        throw InvalidStringError(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('s' | 't') '_'? "
            "Z | 'D',\n where Z matches integers, and ignoring whitespaces.");
        ;
    }
}

void Underlying::direct(sint16 *dir_perm) const {
    for (sint16 i = 0; i < get_parameter().n; i++) {
        dir_perm[permutation_table[i]] = i;
    }
};

Underlying Underlying::left_meet(const Underlying &b) const {
    thread_local sint16 dir_perm_a[MAX_N], dir_perm_b[MAX_N], dir_perm_meet[MAX_N];
    Underlying a_copy = *this;
    Underlying b_copy = b;
    Underlying meet(get_parameter());
    meet.identity();
    a_copy.direct(dir_perm_a);
    b_copy.direct(dir_perm_b);
    meet.direct(dir_perm_meet);
    sint16 n = get_parameter().n, e = get_parameter().e;
    for (sint16 i = 2; i <= n; i++) {
        for (sint16 j = i; j > 2; j--) {
            if (a_copy.is_s_left_divisor(dir_perm_a, j) &&
                b_copy.is_s_left_divisor(dir_perm_b, j)) {
                meet.s_right_multiply(dir_perm_meet, j);
                a_copy.s_left_multiply(dir_perm_a, j);
                b_copy.s_left_multiply(dir_perm_b, j);
            } else {
                goto end_subword;
            }
        }

        if (a_copy.is_t_left_divisor(dir_perm_a, 0) &&
            b_copy.is_t_left_divisor(dir_perm_b, 0)) {

            meet.t_right_multiply(dir_perm_meet, 0);
            a_copy.t_left_multiply(dir_perm_a, 0);
            b_copy.t_left_multiply(dir_perm_b, 0);

            if (a_copy.is_t_left_divisor(dir_perm_a, e - 1) &&
                b_copy.is_t_left_divisor(dir_perm_b, e - 1)) {
                meet.t_right_multiply(dir_perm_meet, e - 1);
                a_copy.t_left_multiply(dir_perm_a, e - 1);
                b_copy.t_left_multiply(dir_perm_b, e - 1);

                for (sint16 j = 3; j < i + 1; j++) {
                    if (a_copy.is_s_left_divisor(dir_perm_a, j) &&
                        b_copy.is_s_left_divisor(dir_perm_b, j)) {
                        meet.s_right_multiply(dir_perm_meet, j);
                        a_copy.s_left_multiply(dir_perm_a, j);
                        b_copy.s_left_multiply(dir_perm_b, j);
                    } else {
                        goto end_subword;
                    }
                }
            }
            goto end_subword;
        }

        for (sint16 k = 1; k < e; k++) {
            if (a_copy.is_t_left_divisor(dir_perm_a, k) &&
                b_copy.is_t_left_divisor(dir_perm_b, k)) {
                meet.t_right_multiply(dir_perm_meet, k);
                a_copy.t_left_multiply(dir_perm_a, k);
                b_copy.t_left_multiply(dir_perm_b, k);
                goto end_subword;
            }
        }
    end_subword:;
    }
    return meet;
}

void Underlying::identity() {
    for (sint16 i = 0; i < get_parameter().n; i++) {
        permutation_table[i] = i;
        coefficient_table[i] = 0;
    }
}

void Underlying::delta() {
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 1; i < n; i++) {
        permutation_table[i] = i;
        coefficient_table[i] = 1;
    }
    permutation_table[0] = 0;
    coefficient_table[0] = Rem(-n + 1, e);
}

bool Underlying::compare(const Underlying &b) const {
    sint16 i;
    for (i = 0; i < get_parameter().n; i++) {
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
    for (i = 0; i < n; i++) {
        f.permutation_table[permutation_table[i]] = i;
        f.coefficient_table[permutation_table[i]] =
            coefficient_table[i] == 0 ? 0 : e - coefficient_table[i];
    }
    return f;
};

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    sint16 i, n = get_parameter().n, e = get_parameter().e;
    for (i = 0; i < n; i++) {
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
    // delta is diagonal, and acts almost homothetically.
    // Therefore conjugating by some power of it does nothing on most
    // coefficients.
    sint16 n = get_parameter().n, e = get_parameter().e;

    // In this case `*this` commutes with delta.
    if (permutation_table[0] == 0) {
        return;
    }

    // Otherwise the two non trivial coefficient are 0 and the i such that
    // `permutation_table[i] == 0`.
    coefficient_table[0] = Rem(coefficient_table[0] + k * n, e);
    for (sint16 i = 1; i < n; i++) {
        if (permutation_table[i] == 0) {
            coefficient_table[i] = Rem(coefficient_table[i] - k * n, e);
            return;
        }
    }
}

void Underlying::randomize() { throw NonRandomizable(); }

std::vector<Underlying> Underlying::atoms() const {
    Parameter p = get_parameter();
    std::vector<Underlying> atoms;
    Underlying atom = Underlying(p);
    sint16 n = p.n, e = p.e;
    for (sint16 i = 2; i <= n - 1; i++) {
        // s_(i+1).
        atom.identity();
        atom.permutation_table[i - 1] = i;
        atom.permutation_table[i] = i - 1;
        atoms.push_back(atom);
    }
    for (sint16 k = 0; k < e; k++) {
        // t_k.
        atom.identity();
        atom.permutation_table[0] = 1;
        atom.permutation_table[1] = 0;
        atom.coefficient_table[1] = k;
        atom.coefficient_table[0] = e - ((k == 0) ? e : k);
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace standard_complex

template <>
IndentedOStream &IndentedOStream::operator<< <standard_complex::EENParameter>(
    const standard_complex::EENParameter &p) {
    p.print(*this);
    return *this;
};

} // namespace cgarside
/**
 * @file dihedral_braid.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for I-series Artin groups (dual Garside structure).
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

#include "dihedral.h"

struct NotBelow {};

namespace garcide::dihedral {

sint16 Underlying::get_parameter() const { return number_of_points; }

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
            throw InvalidStringError("Parameter is too big!\n" + match.str(1) +
                                     " can not be converted to a C++ integer.");
        }
        if (2 <= i) {
            return i;
        } else {
            throw InvalidStringError("Parameter should be at least 2!");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
};

sint16 Underlying::lattice_height() const { return 2; }

Underlying::Underlying(sint16 n)
    : number_of_points(n), type(0), point(0) {}

void Underlying::print(IndentedOStream &os) const {
    if (type == 1) {
        os << "D";
    } else if (type == 2) {
        os << "s" << point;
    }
}

void Underlying::of_string(const std::string &str, size_t &pos) {
    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"(s[\\s\\t]*_?)?[\\s\\t]*(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        point = Rem(std::stoi(match[2]), get_parameter());
        type = 2;
        pos += match[0].length();
    } else {
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('s' '_'?)? Z | 'D',\nwhere Z "
            "matches integers, and ignoring whitespaces."));
    }
}

Underlying Underlying::left_meet(const Underlying &b) const {
    if ((type == 0) || (b.type == 0) ||
        ((type == 2) && (b.type == 2) && (point != b.point))) {
        return Underlying(get_parameter());
    }
    if (type == 1) {
        Underlying c = b;
        return c;
    }
    Underlying c = *this;
    return c;
}

Underlying Underlying::right_meet(const Underlying &b) const {
    return left_meet(b);
}

void Underlying::identity() { type = 0; }

void Underlying::delta() { type = 1; }

bool Underlying::compare(const Underlying &b) const {
    return ((type == b.type) && ((type != 2) || (point == b.point)));
}

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (type == 0) {
        f = b;
    } else if (b.type == 0) {
        f = *this;
    } else if ((type == 2) && (b.type == 2) &&
               ((point - b.point + get_parameter()) % get_parameter() == 1)) {
        f.delta();
    } else {
        throw NotBelow();
    }
    return f;
};

Underlying Underlying::left_complement(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (b.type == 1) {
        if (type == 0) {
            f.delta();
        } else if (type == 2) {
            f.type = 2;
            if (point == get_parameter() - 1) {
                f.point = 0;
            } else {
                f.point = point + 1;
            }
        }
    } else if (b.type == 2) {
        if (type == 0) {
            f = b;
        } else if ((type == 2) && (point == b.point)) {
            f.identity();
        } else {
            throw NotBelow();
        }
    } else {
        if (type != 0) {
            throw NotBelow();
        }
    }
    return f;
};

Underlying Underlying::right_complement(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (b.type == 1) {
        if (type == 0) {
            f.delta();
        } else if (type == 2) {
            f.type = 2;
            if (point == 0) {
                f.point = get_parameter() - 1;
            } else {
                f.point = point - 1;
            }
        }
    } else if (b.type == 2) {
        if (type == 0) {
            f = b;
        } else if ((type == 2) && (point == b.point)) {
            f.identity();
        } else {
            throw NotBelow();
        }
    } else {
        if (type != 0) {
            throw NotBelow();
        }
    }
    return f;
};

Underlying Underlying::delta_conjugate_mut(sint16 k) const {
    Underlying under = Underlying(*this);
    sint16 n = get_parameter();
    if (type != 2) {
        under = *this;
    } else {
        if (k > 0) {
            k = k - k * n;
        }
        under.type = 2;
        under.point = (point - 2 * k) % n;
    }

    return under;
}

void Underlying::randomize() {
    sint16 rand = std::rand() % (get_parameter() + 1);
    if (rand == get_parameter()) {
        type = 0;
    } else if (rand == get_parameter() + 1) {
        type = 1;
    } else {
        type = 2;
        point = rand;
    }
}

std::vector<Underlying> Underlying::atoms() const {
    sint16 n = get_parameter();
    Underlying atom = Underlying(n);
    atom.type = 2;
    std::vector<Underlying> atoms;
    for (sint16 i = 0; i < n; i++) {
        atom.point = i;
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace cgarside::dihedral
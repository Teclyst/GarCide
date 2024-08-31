/**
 * @file dihedral.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for \f$\mathbf I\f$-series Artin groups (dual
 * Garside structure).
 * @version 1.0.0
 * @date 2024-08-31
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

#include "garcide/groups/dihedral.hpp"

namespace garcide::dihedral {

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
        if (2 <= i) {
            return i;
        } else {
            throw InvalidStringError("Parameter should be at least 2!");
        }
    } else {
        throw InvalidStringError(
            std::string("Could not extract an integer from \"str\"!"));
    }
}

Underlying::Underlying(i16 n) : number_of_vertices(n), type(0), vertex(0) {}

void Underlying::print(IndentedOStream &os) const {
    if (type == 1) {
        os << "D";
    } else if (type == 2) {
        os << "s" << vertex;
    }
}

void Underlying::of_string(const std::string &str, size_t &pos) {
    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"(:?s[\\s\\t]*_?)?[\\s\\t]*(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        vertex = rem(std::stoi(match[1]), get_parameter());
        type = 2;
        pos += match[0].length();
    } else {
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('s' '_'?)? Z | 'D',\nwhere Z "
            "matches integers, and ignoring whitespaces."));
    }
}

void Underlying::debug(IndentedOStream &os) const {
    os << "{   ";
    os.Indent(4);
    os << "number_of_vertices:";
    os.Indent(4);
    os << EndLine() << number_of_vertices;
    os.Indent(-4);
    os << EndLine();
    os << "type:";
    os.Indent(4);
    os << EndLine();
    os << type;
    os.Indent(-4);
    os << EndLine();
    os << "vertex:";
    os.Indent(4);
    os << EndLine();
    os << vertex;
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

Underlying Underlying::left_meet(const Underlying &b) const {
    if ((type == 0) || (b.type == 0) ||
        ((type == 2) && (b.type == 2) && (vertex != b.vertex))) {
        return Underlying(get_parameter());
    }
    if (type == 1) {
        Underlying c = b;
        return c;
    }
    Underlying c = *this;
    return c;
}

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (type == 0) {
        f = b;
    } else if (b.type == 0) {
        f = *this;
    } else if ((type == 2) && (b.type == 2) &&
               ((vertex - b.vertex + get_parameter()) % get_parameter() == 1)) {
        f.delta();
    } else {
        throw NotBelow();
    }
    return f;
}

Underlying Underlying::left_complement(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (b.type == 1) {
        if (type == 0) {
            f.delta();
        } else if (type == 2) {
            f.type = 2;
            if (vertex == get_parameter() - 1) {
                f.vertex = 0;
            } else {
                f.vertex = vertex + 1;
            }
        }
    } else if (b.type == 2) {
        if (type == 0) {
            f = b;
        } else if ((type == 2) && (vertex == b.vertex)) {
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
}

Underlying Underlying::right_complement(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (b.type == 1) {
        if (type == 0) {
            f.delta();
        } else if (type == 2) {
            f.type = 2;
            if (vertex == 0) {
                f.vertex = get_parameter() - 1;
            } else {
                f.vertex = vertex - 1;
            }
        }
    } else if (b.type == 2) {
        if (type == 0) {
            f = b;
        } else if ((type == 2) && (vertex == b.vertex)) {
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
}

void Underlying::delta_conjugate_mut(i16 k) {
    if (type == 2) {
        vertex = rem(vertex - 2 * k, get_parameter());
    }
}

void Underlying::randomize() {
    i16 rand = std::rand() % (get_parameter() + 1);
    if (rand == get_parameter()) {
        type = 0;
    } else if (rand == get_parameter() + 1) {
        type = 1;
    } else {
        type = 2;
        vertex = rand;
    }
}

std::vector<Underlying> Underlying::atoms() const {
    i16 n = get_parameter();
    Underlying atom = Underlying(n);
    atom.type = 2;
    std::vector<Underlying> atoms;
    for (i16 i = 0; i < n; i++) {
        atom.vertex = i;
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace garcide::dihedral
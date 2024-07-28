#include "dihedral_braid.h"
#include <iostream>
#include <string>

struct NotBelow {};

namespace cgarside::dihedral {

sint16 Underlying::get_parameter() const { return PresentationParameter; }

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
    : PresentationParameter(n), Type(0), Point(0) {}

void Underlying::print(IndentedOStream &os) const {
    if (Type == 1) {
        os << "D";
    } else if (Type == 2) {
        os << "s" << Point;
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
        Point = Rem(std::stoi(match[2]), get_parameter());
        Type = 2;
        pos += match[0].length();
    } else {
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('s' '_'?)? Z | 'D',\nwhere Z "
            "matches integers, and ignoring whitespaces."));
    }
}

Underlying Underlying::left_meet(const Underlying &b) const {
    if ((Type == 0) || (b.Type == 0) ||
        ((Type == 2) && (b.Type == 2) && (Point != b.Point))) {
        return Underlying(get_parameter());
    }
    if (Type == 1) {
        Underlying c = b;
        return c;
    }
    Underlying c = *this;
    return c;
}

Underlying Underlying::right_meet(const Underlying &b) const {
    return left_meet(b);
}

void Underlying::identity() { Type = 0; }

void Underlying::delta() { Type = 1; }

bool Underlying::compare(const Underlying &b) const {
    return ((Type == b.Type) && ((Type != 2) || (Point == b.Point)));
}

Underlying Underlying::product(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (Type == 0) {
        f = b;
    } else if (b.Type == 0) {
        f = *this;
    } else if ((Type == 2) && (b.Type == 2) &&
               ((Point - b.Point + get_parameter()) % get_parameter() == 1)) {
        f.delta();
    } else {
        throw NotBelow();
    }
    return f;
};

Underlying Underlying::left_complement(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (b.Type == 1) {
        if (Type == 0) {
            f.delta();
        } else if (Type == 2) {
            f.Type = 2;
            if (Point == get_parameter() - 1) {
                f.Point = 0;
            } else {
                f.Point = Point + 1;
            }
        }
    } else if (b.Type == 2) {
        if (Type == 0) {
            f = b;
        } else if ((Type == 2) && (Point == b.Point)) {
            f.identity();
        } else {
            throw NotBelow();
        }
    } else {
        if (Type != 0) {
            throw NotBelow();
        }
    }
    return f;
};

Underlying Underlying::right_complement(const Underlying &b) const {
    Underlying f = Underlying(get_parameter());
    if (b.Type == 1) {
        if (Type == 0) {
            f.delta();
        } else if (Type == 2) {
            f.Type = 2;
            if (Point == 0) {
                f.Point = get_parameter() - 1;
            } else {
                f.Point = Point - 1;
            }
        }
    } else if (b.Type == 2) {
        if (Type == 0) {
            f = b;
        } else if ((Type == 2) && (Point == b.Point)) {
            f.identity();
        } else {
            throw NotBelow();
        }
    } else {
        if (Type != 0) {
            throw NotBelow();
        }
    }
    return f;
};

Underlying Underlying::delta_conjugate_mut(sint16 k) const {
    Underlying under = Underlying(*this);
    sint16 n = get_parameter();
    if (Type != 2) {
        under = *this;
    } else {
        if (k > 0) {
            k = k - k * n;
        }
        under.Type = 2;
        under.Point = (Point - 2 * k) % n;
    }

    return under;
}

void Underlying::randomize() {
    sint16 rand = std::rand() % (get_parameter() + 1);
    if (rand == get_parameter()) {
        Type = 0;
    } else if (rand == get_parameter() + 1) {
        Type = 1;
    } else {
        Type = 2;
        Point = rand;
    }
}

std::vector<Underlying> Underlying::atoms() const {
    sint16 n = get_parameter();
    Underlying atom = Underlying(n);
    atom.Type = 2;
    std::vector<Underlying> atoms;
    for (sint16 i = 0; i < n; i++) {
        atom.Point = i;
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace cgarside::dihedral
#include "dihedral_braid.h"
#include <iostream>
#include <string>

struct NotBelow {};

namespace cgarside::dihedral {

sint16 Underlying::GetParameter() const { return PresentationParameter; }

Underlying::ParameterType
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

sint16 Underlying::LatticeHeight() const { return 2; }

Underlying::Underlying(sint16 n)
    : PresentationParameter(n), Type(0), Point(0) {}

void Underlying::Print(IndentedOStream &os) const {
    if (Type == 1) {
        os << "D";
    } else if (Type == 2) {
        os << "s" << Point;
    }
}

void Underlying::OfString(const std::string &str, size_t &pos) {
    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match, std::regex{"D"},
                          std::regex_constants::match_continuous)) {
        pos += match[0].length();
        Delta();
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"(s[\\s\\t]*_?)?[\\s\\t]*(" +
                                            number_regex + ")"},
                                 std::regex_constants::match_continuous)) {
        Point = Rem(std::stoi(match[2]), GetParameter());
        Type = 2;
        pos += match[0].length();
    } else {
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('s' '_'?)? Z | 'D',\nwhere Z "
            "matches integers, and ignoring whitespaces."));
    }
}

Underlying Underlying::LeftMeet(const Underlying &b) const {
    if ((Type == 0) || (b.Type == 0) ||
        ((Type == 2) && (b.Type == 2) && (Point != b.Point))) {
        return Underlying(GetParameter());
    }
    if (Type == 1) {
        Underlying c = b;
        return c;
    }
    Underlying c = *this;
    return c;
}

Underlying Underlying::RightMeet(const Underlying &b) const {
    return LeftMeet(b);
}

void Underlying::Identity() { Type = 0; }

void Underlying::Delta() { Type = 1; }

bool Underlying::Compare(const Underlying &b) const {
    return ((Type == b.Type) && ((Type != 2) || (Point == b.Point)));
}

Underlying Underlying::Product(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    if (Type == 0) {
        f = b;
    } else if (b.Type == 0) {
        f = *this;
    } else if ((Type == 2) && (b.Type == 2) &&
               ((Point - b.Point + GetParameter()) % GetParameter() == 1)) {
        f.Delta();
    } else {
        throw NotBelow();
    }
    return f;
};

Underlying Underlying::LeftComplement(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    if (b.Type == 1) {
        if (Type == 0) {
            f.Delta();
        } else if (Type == 2) {
            f.Type = 2;
            if (Point == GetParameter() - 1) {
                f.Point = 0;
            } else {
                f.Point = Point + 1;
            }
        }
    } else if (b.Type == 2) {
        if (Type == 0) {
            f = b;
        } else if ((Type == 2) && (Point == b.Point)) {
            f.Identity();
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

Underlying Underlying::RightComplement(const Underlying &b) const {
    Underlying f = Underlying(GetParameter());
    if (b.Type == 1) {
        if (Type == 0) {
            f.Delta();
        } else if (Type == 2) {
            f.Type = 2;
            if (Point == 0) {
                f.Point = GetParameter() - 1;
            } else {
                f.Point = Point - 1;
            }
        }
    } else if (b.Type == 2) {
        if (Type == 0) {
            f = b;
        } else if ((Type == 2) && (Point == b.Point)) {
            f.Identity();
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

Underlying Underlying::DeltaConjugate(sint16 k) const {
    Underlying under = Underlying(*this);
    sint16 n = GetParameter();
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

void Underlying::Randomize() {
    sint16 rand = std::rand() % (GetParameter() + 1);
    if (rand == GetParameter()) {
        Type = 0;
    } else if (rand == GetParameter() + 1) {
        Type = 1;
    } else {
        Type = 2;
        Point = rand;
    }
}

std::vector<Underlying> Underlying::Atoms() const {
    sint16 n = GetParameter();
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
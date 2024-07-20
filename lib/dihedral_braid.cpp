#include "dihedral_braid.h"
#include <iostream>
#include <string>

struct NotBelow {};

namespace CGarside {

sint16 IDualBraidUnderlying::GetParameter() const {
    return PresentationParameter;
}

sint16 IDualBraidUnderlying::LatticeHeight() const { return 2; }

IDualBraidUnderlying::IDualBraidUnderlying(sint16 n)
    : PresentationParameter(n), Type(0), Point(0) {}


void IDualBraidUnderlying::Print(IndentedOStream &os) const {
    if (Type == 1) {
        os << "D";
    } else if (Type == 2) {
        os << "s" << Point;
    }
}

void IDualBraidUnderlying::OfString(const std::string &str, size_t &pos) {
    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"(:?s[\\s\\t]*_)?[\\s\\t]*(" + number_regex + ")"},
                          std::regex_constants::match_continuous)) {
        Point = Rem(std::stoi(match[1]), GetParameter());
        Type = 2;
        pos += match[0].length();
    } else {
        throw InvalidStringError(
            std::string("Could not extract a factor from \"" + str.substr(pos) +
                        "\"!\nA factor should match regex (s _?)? -? [1 - 9] [0 - 9]*. (Ignoring whitespaces)"));
    }
}

IDualBraidUnderlying
IDualBraidUnderlying::LeftMeet(const IDualBraidUnderlying &b) const {
    if ((Type == 0) || (b.Type == 0) ||
        ((Type == 2) && (b.Type == 2) && (Point != b.Point))) {
        return IDualBraidUnderlying(GetParameter());
    }
    if (Type == 1) {
        IDualBraidUnderlying c = b;
        return c;
    }
    IDualBraidUnderlying c = *this;
    return c;
}

IDualBraidUnderlying
IDualBraidUnderlying::RightMeet(const IDualBraidUnderlying &b) const {
    return LeftMeet(b);
}

void IDualBraidUnderlying::Identity() { Type = 0; }

void IDualBraidUnderlying::Delta() { Type = 1; }

bool IDualBraidUnderlying::Compare(const IDualBraidUnderlying &b) const {
    return ((Type == b.Type) && ((Type != 2) || (Point == b.Point)));
}

IDualBraidUnderlying
IDualBraidUnderlying::Product(const IDualBraidUnderlying &b) const {
    IDualBraidUnderlying f = IDualBraidUnderlying(GetParameter());
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

IDualBraidUnderlying
IDualBraidUnderlying::LeftComplement(const IDualBraidUnderlying &b) const {
    IDualBraidUnderlying f = IDualBraidUnderlying(GetParameter());
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

IDualBraidUnderlying
IDualBraidUnderlying::RightComplement(const IDualBraidUnderlying &b) const {
    IDualBraidUnderlying f = IDualBraidUnderlying(GetParameter());
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

IDualBraidUnderlying IDualBraidUnderlying::DeltaConjugate(sint16 k) const {
    IDualBraidUnderlying under = IDualBraidUnderlying(*this);
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

void IDualBraidUnderlying::Randomize() {
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

std::vector<IDualBraidUnderlying> IDualBraidUnderlying::Atoms() const {
    sint16 n = GetParameter();
    IDualBraidUnderlying atom = IDualBraidUnderlying(n);
    atom.Type = 2;
    std::vector<IDualBraidUnderlying> atoms;
    for (sint16 i = 0; i < n; i++) {
        atom.Point = i;
        atoms.push_back(atom);
    }
    return atoms;
}

} // namespace CGarside
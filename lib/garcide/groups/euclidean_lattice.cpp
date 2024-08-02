#include "garcide/groups/euclidean_lattice.hpp"

namespace garcide::euclidean_lattice {

Underlying::Parameter Underlying::parameter_of_string(const std::string &str) {
    std::smatch match;

    if (std::regex_match(str, match,
                         std::regex{"[\\s\\t]*([1-9][0-9]*)[\\s\\t]*"})) {
        // Try to match str with a strictly positive integer.
        try {
            return std::stoull(match[1]);
        } catch (std::out_of_range const &) {
            // Even if str matches a positive integer, it might be too big for
            // C++.
            throw InvalidStringError("Dimension is too big!\n" + match.str(1) +
                                     " cannot be converted to a C++ integer.");
        }
    } else {
        // Otherwise fail and explain why.
        throw InvalidStringError(std::string(
            "Could not extract a strictly positive integer from \"" + str +
            "\"!"));
    }
}

Underlying::Parameter Underlying::get_parameter() const {
    return coordinates.size();
}

Underlying::Underlying(Underlying::Parameter n) : coordinates(n, false) {}

void Underlying::of_string(const std::string &str, size_t &pos) {
    Parameter n = get_parameter();

    std::smatch match;

    if (std::regex_search(str.begin() + pos, str.end(), match,
                          std::regex{"(e_?)?(" + number_regex + ")"},
                          std::regex_constants::match_continuous)) {
        // Try to extract a substring starting at pos that matches the regex.
        sint16 i;
        try {
            i = std::stoi(match[2]);
        } catch (std::out_of_range const &) {
            throw InvalidStringError("Index is too big!\n" + match.str(2) +
                                     " can not be converted to a C++ integer.");
        }
        pos += match[0].length();
        if ((i >= 0) && (i < int(n))) {
            identity();
            coordinates[i] = true;
        } else {
            throw InvalidStringError(
                "Invalid index for canonical base vector!\n" + match.str(2) +
                " is not in [0, " + std::to_string(n) + "[.");
        }
    } else if (std::regex_search(str.begin() + pos, str.end(), match,
                                 std::regex{"D"},
                                 std::regex_constants::match_continuous)) {
        // Else try to match "D".
        pos += match[0].length();
        delta();
    } else {
        // Otherwise fail.
        throw InvalidStringError(std::string(
            "Could not extract a factor from\n\"" + str.substr(pos) +
            "\"!\nA factor should match regex ('e' '_'?)? Z | 'D',\nwhere Z "
            "matches integers and ignoring whitespaces."));
    }
}

void Underlying::debug(IndentedOStream &os) const {
    os << "{   ";
    os.Indent(4);
    os << "coordinates:";
    os.Indent(4);
    os << EndLine();
    os << "[";
    for (size_t i = 0; i < get_parameter() - 1; i++) {
        os << coordinates[i] << ", ";
    }
    os << coordinates[get_parameter() - 1];
    os << "]";
    os.Indent(-8);
    os << EndLine();
    os << "}";
}

void Underlying::print(IndentedOStream &os) const {
    Parameter n = get_parameter();

    Underlying c = Underlying(*this);

    bool is_first = true;

    for (size_t i = 0; i < n; i++) {
        if (at(i)) {
            os << (is_first ? "" : " ") << "e" << i;
            is_first = false;
        }
    }
}

void Underlying::identity() {
    for (size_t i = 0; i < get_parameter(); i++) {
        coordinates[i] = false;
    }
}

void Underlying::delta() {
    for (size_t i = 0; i < get_parameter(); i++) {
        coordinates[i] = true;
    }
}

Underlying Underlying::product(const Underlying &b) const {
    Underlying product(get_parameter());
    for (size_t i = 0; i < get_parameter(); i++) {
        product.coordinates[i] = at(i) ^ b.at(i);
    }
    return product;
}

Underlying Underlying::left_meet(const Underlying &b) const {
    Underlying meet(get_parameter());
    for (size_t i = 0; i < get_parameter(); i++) {
        meet.coordinates[i] = at(i) && b.at(i);
    }
    return meet;
}

std::vector<Underlying> Underlying::atoms() const {
    std::vector<Underlying> atoms;
    Underlying atom(get_parameter());
    for (size_t i = 0; i < get_parameter(); i++) {
        atom.coordinates[i] = true;
        atoms.push_back(atom);
        atom.coordinates[i] = false;
    }
    return atoms;
}

size_t Underlying::hash() const {
    size_t hash = 0;
    for (size_t i = 0; i < get_parameter(); i++) {
        hash <<= 1;
        hash += (at(i) ? 1 : 0);
    }
    return hash;
}

} // namespace garcide::euclidean_lattice
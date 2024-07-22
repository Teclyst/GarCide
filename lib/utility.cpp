#include "utility.h"

namespace cgarside {

EndLine::EndLine(sint16 skip) : lines_to_skip(skip) {}

IndentedOStream::IndentedOStream(std::ostream &os) : indent_level(0), os(os) {}

template <>
IndentedOStream &IndentedOStream::operator<< <EndLine>(const EndLine &el) {
    for (sint16 _ = 0; _ < el.lines_to_skip; _++) {
        os << "\n";
    }
    os << std::endl;
    for (sint16 _ = 0; _ < indent_level; _++) {
        os << " ";
    }
    return *this;
}

} // namespace CGarside
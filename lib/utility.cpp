#include "utility.h"

namespace Utility {

IndentedOStream::IndentedOStream(std::ostream &os) : indent_level(0), os(os) {}

template <>
IndentedOStream &IndentedOStream::operator<< <EndLine>(const EndLine &el) {
    os << std::endl;
    for (sint16 _ = 0; _ < indent_level; _++) {
        os << " ";
    }
    return *this;
}

} // namespace Utility
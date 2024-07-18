#include "utility.h"

namespace Utility {

IndentedOStream::IndentedOStream(std::ostream &os) : indent_level(0), os(os) {}

inline void IndentedOStream::SetIndentLevel(sint16 new_indent_level) {
    indent_level = new_indent_level;
}

inline void IndentedOStream::Indent(sint16 indentation) {
    indent_level += indentation;
}

template <>
IndentedOStream &IndentedOStream::operator<< <EndLine>(const EndLine &el) {
    os << std::endl;
    for (sint16 _ = 0; _ < indent_level; _++) {
        os << " ";
    }
    return *this;
}

} // namespace Utility
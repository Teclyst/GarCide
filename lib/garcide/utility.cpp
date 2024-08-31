/**
 * @file utility.cpp
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Implementation file for a bunch of utility functions used everywhere.
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

#include "garcide/utility.hpp"

namespace garcide {

EndLine::EndLine(i16 skip) : lines_to_skip(skip) {}

IndentedOStream::IndentedOStream(std::ostream &os) : indent_level(0), os(os) {}

template <>
IndentedOStream &IndentedOStream::operator<< <EndLine>(const EndLine &el) {
    for (i16 _ = 0; _ < el.lines_to_skip; _++) {
        os << "\n";
    }
    os << std::endl;
    for (i16 _ = 0; _ < indent_level; _++) {
        os << " ";
    }
    return *this;
}

} // namespace garcide
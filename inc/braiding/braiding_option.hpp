/**
 * @file braiding_option.h
 * @author Matteo Wei (matteo.wei@ens.psl.eu)
 * @brief Header file to configure the braiding executable.
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

#ifndef BRAIDING_OPTION
#define BRAIDING_OPTION

#if BRAIDING_CLASS == 0

#include "garcide/groups/artin.hpp"

#elif BRAIDING_CLASS == 1

#include "garcide/groups/band.hpp"

#elif BRAIDING_CLASS == 2

#include "garcide/groups/octahedral.hpp"

#elif BRAIDING_CLASS == 3

#include "garcide/groups/dihedral.hpp"

#elif BRAIDING_CLASS == 4

#include "garcide/groups/dual_complex.h"

#elif BRAIDING_CLASS == 5

#include "garcide/groups/standard_complex.h"

#elif BRAIDING_CLASS == 6

#include "garcide/groups/euclidean_lattice.hpp"

#endif

namespace braiding {
#if BRAIDING_CLASS == 0

/**
 * @brief The factor class used for _Braiding_.
 * 
 * Depends on preprocessor variable `BRAIDING_CLASS`. Therefore, it isn't
 * necessarily set to `garcide::artin::Factor`.
 */
using Factor = garcide::artin::Factor;

/**
 * @brief The braid class used for _Braiding_.
 * 
 * Depends on preprocessor variable `BRAIDING_CLASS`. Therefore, it isn't
 * necessarily set to `garcide::artin::Braid`.
 */
using Braid = garcide::artin::Braid;

#elif BRAIDING_CLASS == 1

using Factor = garcide::band::Factor;
using Braid = garcide::band::Braid;

#elif BRAIDING_CLASS == 2

using Factor = garcide::octahedral::Factor;
using Braid = garcide::octahedral::Braid;

#elif BRAIDING_CLASS == 3

using Factor = garcide::dihedral::Factor;
using Braid = garcide::dihedral::Braid;

#elif BRAIDING_CLASS == 4

using Factor = garcide::dual_complex::Factor;
using Braid = garcide::dual_complex::Braid;

#elif BRAIDING_CLASS == 5

using Factor = garcide::standard_complex::Factor;
using Braid = garcide::standard_complex::Braid;

#elif BRAIDING_CLASS == 6

using Factor = garcide::euclidean_lattice::Factor;
using Braid = garcide::euclidean_lattice::Braid;

#endif

}

#endif
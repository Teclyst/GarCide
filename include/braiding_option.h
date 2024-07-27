#ifndef BRAIDING_OPTION
#define BRAIDING_OPTION

#if BRAIDING_CLASS == 0

#include "artin_braid.h"

#elif BRAIDING_CLASS == 1

#include "band_braid.h"

#elif BRAIDING_CLASS == 2

#include "octahedral_braid.h"

#elif BRAIDING_CLASS == 3

#include "dihedral_braid.h"

#elif BRAIDING_CLASS == 4

#include "dual_complex_reflection.h"

#elif BRAIDING_CLASS == 5

#include "standard_complex_reflection.h"

#endif

namespace braiding {
#if BRAIDING_CLASS == 0

using Factor = cgarside::artin::Factor;
using Braid = cgarside::artin::Braid;

#elif BRAIDING_CLASS == 1

using Factor = cgarside::band::Factor;
using Braid = cgarside::band::Braid;


#elif BRAIDING_CLASS == 2

using Factor = cgarside::octahedral::Factor;
using Braid = cgarside::octahedral::Braid;

#elif BRAIDING_CLASS == 3

using Factor = cgarside::dihedral::Factor;
using Braid = cgarside::dihedral::Braid;

#elif BRAIDING_CLASS == 4

using Factor = cgarside::dual_complex::Factor;
using Braid = cgarside::dual_complex::Braid;

#elif BRAIDING_CLASS == 5

using Factor = cgarside::standard_complex::Factor;
using Braid = cgarside::standard_complex::Braid;


#endif

}

#endif
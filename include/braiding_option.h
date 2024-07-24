#ifndef BRAIDING_OPTION
#define BRAIDING_OPTION

#if BRAIDING_CLASS == 0

#include "artin_braid.h"

#endif

#if BRAIDING_CLASS == 1

#include "band_braid.h"

#endif

namespace braiding {
#if BRAIDING_CLASS == 0

using Factor = cgarside::artin::ArtinFactor;
using Braid = cgarside::artin::ArtinBraid;

#endif

#if BRAIDING_CLASS == 1

using Factor = cgarside::band::BandBraidFactor;
using Braid = cgarside::band::BandBraid;

#endif

}

#endif
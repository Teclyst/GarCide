set(HEADERS_PATH "${GarCide_SOURCE_DIR}/inc/garcide/")
set(HEADERS_LIST utility.h 
    groups/artin.h 
    groups/band.h 
    groups/octahedral.h 
    groups/dihedral.h 
    groups/dual_complex.h 
    groups/standard_complex.h 
    )
list(TRANSFORM HEADERS_LIST PREPEND ${HEADERS_PATH})

add_library(
    garcide
    garcide/utility.cpp
    garcide/groups/artin.cpp
    garcide/groups/band.cpp
    garcide/groups/octahedral.cpp
    garcide/groups/dihedral.cpp
    garcide/groups/dual_complex.cpp
    garcide/groups/standard_complex.cpp
    ${HEADERS_LIST})
target_include_directories(garcide PRIVATE ../inc)

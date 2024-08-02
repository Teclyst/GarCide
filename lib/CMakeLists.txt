set(HEADERS_PATH "${GarCide_SOURCE_DIR}/inc/garcide/")
set(HEADERS_LIST
    utility.h 
    groups/artin.h 
    groups/band.h 
    groups/octahedral.h 
    groups/dihedral.h 
    groups/dual_complex.h 
    groups/standard_complex.h
    groups/euclidean_lattice.hpp
)
list(TRANSFORM HEADERS_LIST PREPEND ${HEADERS_PATH})

set(RANDOMIZE_AS_WORD FALSE CACHE BOOL "If enabled, random braids are produced by taking random words in the atoms.")

add_library(
    garcide
    garcide/utility.cpp
    garcide/groups/artin.cpp
    garcide/groups/band.cpp
    garcide/groups/octahedral.cpp
    garcide/groups/dihedral.cpp
    garcide/groups/dual_complex.cpp
    garcide/groups/standard_complex.cpp
    garcide/groups/euclidean_lattice.cpp
    ${HEADERS_LIST}
)
target_include_directories(garcide PRIVATE ../inc)
target_compile_options(garcide PRIVATE -Wall -Wextra -Wpedantic)

if (${USE_PAR} AND ${TBB_FOUND})
    target_compile_definitions(garcide PRIVATE -DUSE_PAR)
    target_link_libraries(garcide PRIVATE TBB::tbb)
endif()

if (${RANDOMIZE_AS_WORD})
    target_compile_definitions(garcide PRIVATE -DRANDOMIZE_AS_WORD)
endif()
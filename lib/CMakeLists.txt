# Define header files.
set(HEADERS_PATH "${GarCide_SOURCE_DIR}/inc/garcide/")
set(HEADERS_LIST
    utility.hpp
    groups/artin.hpp 
    groups/band.hpp 
    groups/octahedral.hpp 
    groups/dihedral.hpp 
    groups/dual_complex.hpp 
    groups/standard_complex.hpp
    groups/euclidean_lattice.hpp
)
list(TRANSFORM HEADERS_LIST PREPEND ${HEADERS_PATH})

# Configure the garcide library target.
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

# Link TBB if it is present and desired.
if (${USE_PAR} AND ${TBB_FOUND})
    target_compile_definitions(garcide PRIVATE -DUSE_PAR)
    target_link_libraries(garcide PRIVATE TBB::tbb)
endif()

# Define the RANDOMIZE_AS_WORD preprocessor variable if asked to.
if (${RANDOMIZE_AS_WORD})
    target_compile_definitions(garcide PRIVATE -DRANDOMIZE_AS_WORD)
endif()
set(HEADERS_PATH "${GarCide_SOURCE_DIR}/inc/")
set(HEADERS_LIST
    braiding/braiding.h
    garcide/centralizer.h
    garcide/sliding_circuits.h
)
list(TRANSFORM HEADERS_LIST PREPEND ${HEADERS_PATH})

set(USE_FOR_BRAIDING "ARTIN" CACHE STRING "Sets the group that is used for Braiding.")

add_executable(
    braiding.exe
    braiding.cpp
    braiding_main.cpp
    ${HEADERS_LIST}
)
target_include_directories(braiding.exe PRIVATE ../inc)
target_link_libraries(braiding.exe PRIVATE garcide)
target_compile_options(braiding.exe PRIVATE -Wall -Wextra -Wpedantic)

if (${USE_FOR_BRAIDING} STREQUAL "ARTIN")
    message("Braiding will use standard braids (classic Garside structure)!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=0)
elseif (${USE_FOR_BRAIDING} STREQUAL "BAND")
    message("Braiding will use standard braids (dual Garside structure)!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=1)
elseif (${USE_FOR_BRAIDING} STREQUAL "OCTAHEDRAL")
    message("Braiding will use B-series Artin groups (dual Garside structure)!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=2)
elseif (${USE_FOR_BRAIDING} STREQUAL "DIHEDRAL")
    message("Braiding will use I-series Artin groups (dual Garside structure)!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=3)
elseif (${USE_FOR_BRAIDING} STREQUAL "DUAL_COMPLEX")
    message("Braiding will use B(e, e, n) (dual Garside structure)!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=4)
elseif (${USE_FOR_BRAIDING} STREQUAL "STANDARD_COMPLEX")
    message("Braiding will use B(e, e, n) (semi-classic Garside structure)!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=5)
elseif (${USE_FOR_BRAIDING} STREQUAL "EUCLIDEAN_LATTICE")
    message("Braiding will use Z ^ n!")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=6)
else()
    message("${USE_FOR_BRAIDING} is not an existing option!\nBraiding will use standard braids (classic Garside structure).")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=0)
endif()

if (${USE_PAR} AND ${TBB_FOUND})
    target_compile_definitions(braiding.exe PRIVATE -DUSE_PAR)
    target_link_libraries(braiding.exe PRIVATE TBB::tbb)
endif()
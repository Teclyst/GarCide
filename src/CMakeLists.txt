# Define header files.
set(HEADERS_PATH "${GarCide_SOURCE_DIR}/inc/")
set(HEADERS_LIST
    braiding/braiding.hpp
    garcide/centralizer.hpp
    garcide/sliding_circuits.hpp
)
list(TRANSFORM HEADERS_LIST PREPEND ${HEADERS_PATH})

# Define the default value (i.e. at the first compilation) of USE_FOR_BRAIDING.
set(USE_FOR_BRAIDING "ARTIN" CACHE STRING "Sets the group that is used for Braiding.")

# Configure the braiding.exe executable target.
add_executable(
    braiding.exe
    braiding.cpp
    braiding_main.cpp
    ${HEADERS_LIST}
)
target_include_directories(braiding.exe PRIVATE ../inc)
target_link_libraries(braiding.exe PRIVATE garcide)
target_compile_options(braiding.exe PRIVATE -Wall -Wextra -Wpedantic)

# Choose which class to use according to the value of USE_FOR_BRAIDING.
message("${BOLD_MAGENTA}Find the group to use for Braiding${COLOUR_RESET}")
if (${USE_FOR_BRAIDING} STREQUAL "ARTIN")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=0)
elseif (${USE_FOR_BRAIDING} STREQUAL "BAND")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=1)
elseif (${USE_FOR_BRAIDING} STREQUAL "OCTAHEDRAL")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=2)
elseif (${USE_FOR_BRAIDING} STREQUAL "DIHEDRAL")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=3)
elseif (${USE_FOR_BRAIDING} STREQUAL "DUAL_COMPLEX")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=4)
elseif (${USE_FOR_BRAIDING} STREQUAL "STANDARD_COMPLEX")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=5)
elseif (${USE_FOR_BRAIDING} STREQUAL "EUCLIDEAN_LATTICE")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=6)
else()
    message(FATAL_ERROR "Invalid option ${BOLD_RED}${USE_FOR_BRAIDING}${COLOUR_RESET}${RED} for ${BOLD_RED}USE_FOR_BRAIDING${COLOUR_RESET}")
    target_compile_definitions(braiding.exe PRIVATE -DBRAIDING_CLASS=0)
endif()
message("-- ${BOLD}USE_FOR_BRAIDING${COLOUR_RESET} is set to ${BOLD}${USE_FOR_BRAIDING}${COLOUR_RESET}")

# Link TBB if it is present and desired.
if (${USE_PAR} AND ${TBB_FOUND})
    target_compile_definitions(braiding.exe PRIVATE -DUSE_PAR)
    target_link_libraries(braiding.exe PRIVATE TBB::tbb)
endif()

# Print this before default CMake output messages.
message("${BOLD_MAGENTA}Generate build files${COLOUR_RESET}")
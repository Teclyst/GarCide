# Minimal CMake version needed.
cmake_minimum_required(VERSION 3.22)

# Declare standards. C++ 17+ is needed for parallel standard library algorithms.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Project declaration.
project(GarCide VERSION 0.1)

# Dependencies.
find_package(TBB)
find_package(Doxygen)

# Default cache options.
set(USE_PAR TRUE CACHE BOOL "Enable parallelism.")
set(GENERATE_DOC TRUE CACHE BOOL "Generate documentation.")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose between release and debug.")
set(DOXYGEN_QUIET NO CACHE STRING "Ask Doxygen to be quiet.")
set(DOXYGEN_WARNINGS NO CACHE STRING "Disable Doxygen warnings.")
set(RANDOMIZE_AS_WORD FALSE CACHE BOOL "If enabled, random braids are produced by taking random words in the atoms.")

# Colours and formating.
string(ASCII 27 ESC)
set(COLOUR_RESET "${ESC}[m")
set(BOLD "${ESC}[1m")
set(RED "${ESC}[31m")
set(GREEN "${ESC}[32m")
set(BOLD_RED "${ESC}[1;31m")
set(BOLD_GREEN "${ESC}[1;32m")
set(BOLD_MAGENTA "${ESC}[1;35m")

# Doxygen configuration.
set(DOXYGEN_HTML_EXTRA_STYLESHEET "doxygen-awesome-css/doxygen-awesome.css")
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_DISABLE_INDEX NO)
set(DOXYGEN_FULL_SIDEBAR NO)
set(DOXYGEN_HTML_COLORSTYLE LIGHT)
set(DOXYGEN_HTML_HEADER "doc/header.html")
set(DOXYGEN_HTML_EXTRA_FILES "doxygen-awesome-css/doxygen-awesome-darkmode-toggle.js" "doxygen-awesome-css/doxygen-awesome-interactive-toc.js")
set(DOXYGEN_SORT_MEMBER_DOCS YES)
set(DOXYGEN_SORT_BRIEF_DOCS YES)
set(DOXYGEN_USE_MATHJAX YES)

# Display which dependencies were found and compile documentation.
message("${BOLD_MAGENTA}Look for dependencies${COLOUR_RESET}")
if (${USE_PAR})
    if (${TBB_FOUND})
        message("${GREEN}-- Found ${BOLD_GREEN}TBB${COLOUR_RESET}")
    else()
        message("${RED}-- Could not find ${BOLD_RED}TBB${COLOUR_RESET},${RED}, parallelism is disabled")
    endif()
endif()
if (${GENERATE_DOC})
    if (${DOXYGEN_FOUND})
        doxygen_add_docs(
            doxygen
            ${PROJECT_SOURCE_DIR}/inc
            ${PROJECT_SOURCE_DIR}/lib
            ${PROJECT_SOURCE_DIR}/src
            ${PROJECT_SOURCE_DIR}/doc/main_page.dox
            ALL
        )
        message("${GREEN}-- Found ${BOLD_GREEN}Doxygen${COLOUR_RESET}")
    else()
        message("${RED}-- Could not find ${BOLD_RED}Doxygen${RED}${COLOUR_RESET}, documentation generation is disabled${COLOUR_RESET}")
    endif()
endif()

# Specify where the other CMakeLists.txt files are.
add_subdirectory(lib)
add_subdirectory(src)
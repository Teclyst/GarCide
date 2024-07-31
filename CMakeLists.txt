cmake_minimum_required(VERSION 3.22)
project(GarCide)

find_package(TBB)

set(USE_PAR TRUE CACHE BOOL "Enable parallelism.")

if (${USE_PAR})
    if (${TBB_FOUND})
        message("Found TBB!")
    else()
        message("Could not find TBB! Parallelism is disabled.")
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(lib)

add_subdirectory(src)
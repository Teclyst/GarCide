set(HEADERS_PATH "${GarCide_SOURCE_DIR}/inc/")
set(HEADERS_LIST
    braiding/braiding.h
    garcide/centralizer.h
    garcide/sliding_circuits.h
)
list(TRANSFORM HEADERS_LIST PREPEND ${HEADERS_PATH})

add_executable(
    braiding.exe
    braiding.cpp
    braiding_main.cpp
    ${HEADERS_LIST}
)
target_include_directories(braiding.exe PRIVATE ../inc)
target_link_libraries(braiding.exe PRIVATE garcide)
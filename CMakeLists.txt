cmake_minimum_required(VERSION 3.28)

project(MatlibMini LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    message(STATUS "ccache found: ${CCACHE_PROGRAM}")
else()
    message(WARNING "ccache not found — consider installing it")
endif()

set(CMAKE_BUILD_TYPE Release CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)

include(FetchContent)

find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found locally, fetching from GitHub...")
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.5.0
    )
    FetchContent_MakeAvailable(Catch2)
endif()

find_package(TBB QUIET)
if(NOT TBB_FOUND)
    message(STATUS "TBB not found locally, fetching from GitHub...")
    FetchContent_Declare(
        oneTBB
        GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
        GIT_TAG        v2022.3.0
        GIT_SHALLOW    ON
    )
    set(TBB_TEST OFF CACHE INTERNAL "")
    FetchContent_MakeAvailable(oneTBB)
endif()

add_executable(main src/main.cc)
target_include_directories(main PRIVATE include)
target_link_libraries(main PRIVATE TBB::tbb)

add_executable(test_matrix tests/test_matrix.cc)
target_include_directories(test_matrix PRIVATE include)
target_link_libraries(test_matrix PRIVATE Catch2::Catch2WithMain TBB::tbb)

include(CTest)
add_test(NAME MatrixTests COMMAND test_matrix)

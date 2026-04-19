cmake_minimum_required(VERSION 3.28)

project(MatlibMini LANGUAGES CXX)

add_executable(tests src/core/test.cc)
target_include_directories(tests PRIVATE include)

include(CTest)
add_test(NAME tests COMMAND tests)

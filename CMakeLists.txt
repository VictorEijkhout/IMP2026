cmake_minimum_required(VERSION 3.14)
project(decomposition LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library
add_library(decomp INTERFACE)
target_include_directories(decomp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Unit tests (Catch2 v3, found via CMAKE_PREFIX_PATH, e.g. `module load catch2/3.16`)
find_package(Catch2 3 REQUIRED)
enable_testing()
add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE decomp Catch2::Catch2WithMain)
target_compile_options(tests PRIVATE -Wall -Wextra)

include(Catch)
catch_discover_tests(tests)

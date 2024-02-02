cmake_minimum_required(VERSION 3.28.1)
project(linalg 
        VERSION 1.0
        DESCRIPTION "A Linear Algebra Library"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)


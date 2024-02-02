# cmake_minimum_required(VERSION 3.28.2)
# project(linalg_test VERSION 1.0)

include_directories(${PROJECT_SOURCE_DIR}/include/linalg)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_executable(
  unit_tests
  unit_tests.cpp
)


target_link_libraries(
  unit_tests 
  GTest::gtest_main
  linalg-lib 
)

include(GoogleTest)
gtest_discover_tests(unit_tests)

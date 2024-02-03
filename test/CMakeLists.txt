include_directories(${PROJECT_SOURCE_DIR}/include/linalg)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  unit_tests
  AllTests.cpp
  MatrixTest.cpp
  EquationSolverTest.cpp
)


target_link_libraries(
  unit_tests 
  GTest::gtest_main
  linalg-lib 
)

include(GoogleTest)
gtest_discover_tests(unit_tests)

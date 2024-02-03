include_directories(${PROJECT_SOURCE_DIR}/include/linalg)
# file(GLOB SOURCES "linalg/*.cpp")
set(SOURCES linalg/Matrix.cpp linalg/Solution.cpp linalg/Vectors.cpp)
add_library(linalg-lib SHARED ${SOURCES})

add_executable (demo linalg/Demo.cpp)
add_executable (theorems linalg/Theorems.cpp)

target_link_libraries (demo linalg-lib)
target_link_libraries(theorems linalg-lib)

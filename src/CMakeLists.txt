include_directories(${PROJECT_SOURCE_DIR}/include/linalg)
file(GLOB SOURCES "linalg/*.cpp")
add_library(linalg-lib SHARED ${SOURCES})

add_executable (demo linalg/Main.cpp)
target_link_libraries (demo linalg-lib)

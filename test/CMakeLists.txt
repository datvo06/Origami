
cmake_minimum_required(VERSION 3.10)

# Project name
project(MemoryBasedGraphTest)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(. ../src/common ../src/graph ../src/StringTokenizer)

# Source files
set(COMMON_HEADERS ../src/common/*.h)
set(GRAPH_HEADERS ../src/graph/*.h)
set(TOKEN_HEADERS ../src/StringTokenizer/*.h)
set(SRC_FILES graph_test.cpp ../src/StringTokenizer/StringTokenizer.cpp)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -Wall -L/usr/local/lib")

# Add executable
add_executable(graph_test ${SRC_FILES})

# Link libraries (if needed)
# target_link_libraries(graph_test /path/to/required/libraries)

# Custom target for all
add_custom_target(all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target graph_test
    DEPENDS graph_test
)

# Custom clean target
add_custom_target(clean
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/graph_test
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/*.o
)

# Build rules for the StringTokenizer library
add_subdirectory(../src/StringTokenizer ${CMAKE_BINARY_DIR}/StringTokenizer)

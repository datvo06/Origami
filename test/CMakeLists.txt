
cmake_minimum_required(VERSION 3.10)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project name
project(MemoryBasedGraphTest)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
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
# Build rules for the StringTokenizer library
add_subdirectory(../src/StringTokenizer ${CMAKE_BINARY_DIR}/StringTokenizer)


cmake_minimum_required(VERSION 3.10)

# Project name
project(StringTokenizer)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable the generation of compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -ansi -Wall")

# Source files
set(SOURCE_FILES StringTokenizer.cpp strtoktest.cpp)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add library for StringTokenizer
add_library(StringTokenizer STATIC StringTokenizer.cpp)

# Add executable for StrTokTest
add_executable(strtoktest strtoktest.cpp)

# Link the StringTokenizer library to StrTokTest executable
target_link_libraries(strtoktest StringTokenizer)


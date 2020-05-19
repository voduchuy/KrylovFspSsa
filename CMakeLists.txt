CMAKE_MINIMUM_REQUIRED(VERSION 2.7)

PROJECT(KrylovFSPSSA LANGUAGES Fortran)
enable_language(Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_Fortran_FLAGS "-fpic -lblas -llapack")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Set the folders to include
include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)

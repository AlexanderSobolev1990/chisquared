cmake_minimum_required(VERSION 3.7)
project(chi_squared LANGUAGES Fortran)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

add_compile_options(
    -Wall
    "$<$<CONFIG:RELEASE>:-O3>"
    "$<$<CONFIG:DEBUG>:-O0;-g>"
)

set(SOURCES
    src/CHI_SQUARED.F
    src/TEST_CHI_SQUARED.F
)

add_executable(${PROJECT_NAME} ${SOURCES})

cmake_minimum_required(VERSION 3.10)
project(FFTWithTests)

set(CMAKE_CXX_STANDARD 20)

add_library(fft STATIC fft.cpp)

option(BUILD_TESTS "Build unit tests Google Test" OFF)

if(BUILD_TESTS)
    message(STATUS "Building with Google Test")

    find_library(GTEST_LIB gtest REQUIRED)
    find_library(GTEST_MAIN_LIB gtest_main REQUIRED)

    enable_testing()
    add_executable(test_fft test/test_fft.cpp)
    target_link_libraries(test_fft fft ${GTEST_LIB} ${GTEST_MAIN_LIB} pthread)

    add_test(NAME FFTTest COMMAND test_fft)

else()
    message(STATUS "Building main program without tests")

    add_executable(fft_main main.cpp)
    target_link_libraries(fft_main fft)
endif()

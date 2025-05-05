cmake_minimum_required(VERSION 3.10)
project(fft)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


add_executable(fft
    main.cpp
    fft.cpp
    fft.h
)
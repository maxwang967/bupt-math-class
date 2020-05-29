cmake_minimum_required(VERSION 3.15)
project(project_fft)

set(CMAKE_CXX_STANDARD 14)

add_executable(project_fft main.cpp Complex.cpp Complex.h my_fft.cpp my_fft.h)
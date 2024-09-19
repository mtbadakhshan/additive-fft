cmake_minimum_required(VERSION 3.10)
project(additive-fft
    DESCRIPTION "Library for additive FFT")

add_executable(benchmark benchmark.cpp)

target_include_directories(
    benchmark
    
    PUBLIC
    depends/libff
)
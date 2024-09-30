cmake_minimum_required(VERSION 3.10)
project(additive-fft
    DESCRIPTION "Library for additive FFT")

# Enable C++11 or higher
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Find the pkg-config package
find_package(PkgConfig REQUIRED)

# Use pkg-config to find libsodium
pkg_check_modules(LIBSODIUM REQUIRED libsodium)

# Add Libff 
add_subdirectory(depends/libff)

add_executable(main main.cpp)

target_include_directories(
    main

    PRIVATE
    depends/libff
    ${LIBSODIUM_INCLUDE_DIRS}
    .
)

target_link_libraries(
    main 
    PRIVATE 
    
    ${LIBSODIUM_LIBRARIES}
    ff
)

# Add Google Benchmark
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable GTest tests in Benchmark" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable testing of the Benchmark" FORCE)

add_subdirectory(depends/benchmark) # Assuming benchmark is cloned in a subdirectory named 'benchmark'

add_executable(run_benchmark benchmark.cpp)


target_include_directories(
    run_benchmark

    PRIVATE
    depends/libff
    depends/benchmark/include
    ${LIBSODIUM_INCLUDE_DIRS}
    .
)

target_link_libraries(
    run_benchmark 
    PRIVATE 
    
    ${LIBSODIUM_LIBRARIES}
    ff
    benchmark::benchmark
)


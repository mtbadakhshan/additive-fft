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


add_subdirectory(depends/libff)

add_executable(benchmark benchmark.cpp)

target_include_directories(
    benchmark

    PRIVATE
    depends/libff
    ${LIBSODIUM_INCLUDE_DIRS}
    .
)


target_link_libraries(
    benchmark 
    PRIVATE 
    
    ${LIBSODIUM_LIBRARIES}
    ff
)




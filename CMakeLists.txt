project(Identification)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE debug)

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

find_package(ROOT REQUIRED)

add_executable(identify src/main.cc)
target_link_libraries(identify PRIVATE ROOT::Tree ROOT::RIO ROOT::MathMore ROOT::ProofBench ROOT::Postscript )
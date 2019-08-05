cmake_minimum_required(VERSION 3.12.0)

project(standalone_operators LANGUAGES CXX)

add_executable(standalone main.cc)
target_compile_features(standalone PRIVATE cxx_std_17)

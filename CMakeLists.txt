cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0148 NEW)
project(radiation)

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(radiation src/module.cpp)
install(TARGETS radiation DESTINATION .)

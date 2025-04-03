cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0148 NEW)
project(em_radiation)

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(em_radiation src/c_api.cpp)
install(TARGETS em_radiation DESTINATION .)

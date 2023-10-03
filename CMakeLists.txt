set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS -pthread)

include_directories(src)

cmake_minimum_required(VERSION 3.4...3.18)
project(_core LANGUAGES CXX)

add_subdirectory(extern/pybind11)
pybind11_add_module(_core
        src/python_bindings.cpp
        src/OctoAssociator.cpp
        src/OctoAssociator.h
        src/types.h
        src/VelocityModel.cpp
        src/VelocityModel.h
        src/types.cpp)

install(TARGETS _core DESTINATION pyocto)

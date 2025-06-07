cmake_minimum_required(VERSION 3.29)
project(OOP_SPICE)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(OOP_SPICE
        circuit_phase1.h
    main.cpp
        circuit_input.cpp
        mna_builder.cpp
        linear_solver.cpp
        linear_solver.h)

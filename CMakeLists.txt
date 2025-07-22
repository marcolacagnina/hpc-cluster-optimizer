cmake_minimum_required(VERSION 3.20)
project(LJOptimizer)

set(CMAKE_CXX_STANDARD 20)


find_package(OpenMP REQUIRED)

add_executable(LJOptimizer 
    src/main.cpp
    include/point3d.hpp
    include/point3d.cpp
    include/cluster.hpp
    include/cluster.cpp
    src/energy.hpp
    src/differential_evolution.hpp
    src/differential_evolution.cpp
    src/runner.hpp
    src/runner.cpp
)

target_link_libraries(LJOptimizer PRIVATE OpenMP::OpenMP_CXX)
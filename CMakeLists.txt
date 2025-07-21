cmake_minimum_required(VERSION 3.26)
project(LJOptimizer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")

link_directories("/opt/homebrew/opt/libomp/lib")

include_directories("/opt/homebrew/opt/libomp/include")

set(CMAKE_CXX_STANDARD 20)


add_executable(LJOptimizer 
	src/main.cpp
        include/point3d.hpp
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



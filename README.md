# An HPC implementation of a revisited DE global optimization for LJ Cluster
This project is a solver designed to optimize the Geometry of the Lennard-Jones atoms clusters. 

## Build with CMake
This project uses CMake as its build system, follow the steps below to build the project using it.

### Prerequisites

- CMake (version 3.20 or higher)
- OpenMP

### Building the project
To build the executable, make sure your default compiler has access to OpenMP.
In the following commands you can see an example of using a specific compiler g++15 for MacOS.

Then run:
```bash
$ mkdir build
$ cd build
$ CC=gcc-15 CXX=g++-15 cmake .. -D CMAKE_BUILD_TYPE=Release
$ make
```

The executable is created into the folder `build`, it is called `LJOptimizer`.

### Run the executables
Make sure you specify the right number of OpenMP threads `OMP_NUM_THREADS` according to your machine

`LJOptimizer` can be executed through
```bash
$ OMP_NUM_THREADS=8 ./LJOptimizer
```

## Output
The output is a `.xyz` file corresponding to the best configuration of atoms found.
It is suggested to view the output file via VMD or PyMOL, which is created after execution within the `build` folder.


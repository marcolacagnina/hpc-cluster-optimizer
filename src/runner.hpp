#ifndef LJOPTIMIZER_RUNNER_HPP
#define LJOPTIMIZER_RUNNER_HPP

#include <vector>
#include <string>
#include <functional>
#include "../include/point3d.hpp"
#include "differential_evolution.hpp"


// Utility function to save an .xyz file directly from an atom array.
void saveAtomsToXYZ(const std::vector<Point3D>& atoms, const std::string& filename, const std::string& comment = "");


// Utility Function to run <k_runs> times the algorithm and select the BEST result.
OptimizationResult runMultipleDE(
        int k_runs,
        int NP,
        int generations,
        int N_atoms,
        double D,
        const std::function<double(const std::vector<Point3D>&)>& energyFunc
);


#endif //LJOPTIMIZER_RUNNER_HPP

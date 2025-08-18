#ifndef DIFFERENTIAL_EVOLUTION_HPP
#define DIFFERENTIAL_EVOLUTION_HPP

#include "../include/cluster.hpp"
#include <vector>
#include <functional>
#include "energy.hpp"


struct OptimizationResult {
    std::vector<Point3D> bestConfiguration;
    double bestEnergy;
};



OptimizationResult differentialEvolution_Adaptive(
        int NP,                        
        int generations,               // Number of generations
        int N_atoms,                   // N
        double D,                     
        const std::function<double(const std::vector<Point3D>&)>& energyFunc,
        double CR = 0.90               // Crossover rate
);


#endif // DIFFERENTIAL_EVOLUTION_HPP


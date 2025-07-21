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


// Algoritmo di Differential Evolution di GEMINI
OptimizationResult differentialEvolution_Adaptive(
        int NP,                        // Numero individui
        int generations,               // Numero generazioni
        int N_atoms,                   // Numero atomi nel cluster
        double D,                      // Lato della griglia iniziale
        const std::function<double(const std::vector<Point3D>&)>& energyFunc,
        double CR = 0.90               // Crossover rate
);


#endif // DIFFERENTIAL_EVOLUTION_HPP


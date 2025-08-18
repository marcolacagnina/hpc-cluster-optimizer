#include "runner.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>



void saveAtomsToXYZ(const std::vector<Point3D>& atoms, const std::string& filename, const std::string& comment) {
    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    outFile << atoms.size() << "\n";
    outFile << comment << "\n";
    for (const auto& atom : atoms) {
        outFile << std::fixed << std::setprecision(5)
                << "X " << atom.x << " " << atom.y << " " << atom.z << "\n";
    }
}


OptimizationResult runMultipleDE(
        int k_runs,
        int NP,
        int generations,
        int N_atoms,
        double D,
        const std::function<double(const std::vector<Point3D>&)>& energyFunc
) {
    OptimizationResult bestResult;
    bestResult.bestEnergy = std::numeric_limits<double>::infinity();

    std::cout << "Starting " << k_runs << " DE runs..." << std::endl;

    for (int i = 0; i < k_runs; ++i) {

        OptimizationResult currentResult = differentialEvolution_Adaptive(
                NP, generations, N_atoms, D, energyFunc
        );

        if (currentResult.bestEnergy < bestResult.bestEnergy) {
            bestResult = currentResult;
            std::cout << "  -> New best found in run " << (i + 1)
                      << "! Energy: " << bestResult.bestEnergy << std::endl;
        }
    }

    std::cout << "\nAll runs completed. Best energy found: " << bestResult.bestEnergy<< std::endl;

    if (!bestResult.bestConfiguration.empty()) {
        std::string comment = "Best configuration from " + std::to_string(k_runs) +
                              " DE runs. Energy: " + std::to_string(bestResult.bestEnergy);
        saveAtomsToXYZ(bestResult.bestConfiguration, "best_cluster.xyz", comment);
        std::cout << "Best configuration saved to best_cluster.xyz" << std::endl;
    }

    return bestResult;
}

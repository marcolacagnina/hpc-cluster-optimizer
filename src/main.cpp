#include <iomanip>
#include "energy.hpp"
#include "runner.hpp"
#include "differential_evolution.hpp"


int main(int argc, char** argv) {

    int N = 22;
    double D = 1.2;
    int generations = 200000;
    int k_runs = 10;
    int NP = 250;

    std::cout << "Running DE with " << k_runs << " runs for N = " << N << std::endl;

    OptimizationResult result = runMultipleDE(k_runs, NP, generations, N, D, Energy::computeLJEnergy_SIMD);

    std::cout << "Best energy found: " << result.bestEnergy << std::endl;

    return 0;
}




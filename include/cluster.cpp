#include "cluster.hpp"
#include <iomanip>
#include <cmath>


// Constructor
Cluster::Cluster(int N, double D, const std::function<double(const std::vector<Point3D>&)>& energyFunc, std::mt19937& rng)
        : energyFunction(energyFunc)
{
    // Use the provided generator, avoiding the creation of a new std::random_device.
    generateRandomBoxCluster(N, D, rng);
    Energy = energyFunction(atoms);
}


// Cluster Initialization Function.
// Use a random generator, provided by the calling thread, for performance optimization reasons.
void Cluster::generateRandomBoxCluster(int N, double D, std::mt19937& rng) {
    double noise_scale = 0.1;
    atoms.clear(); // Ensures the vector is empty before filling it
    atoms.reserve(N);

    // Compute the cube root of N.
    // Goal: Compute the length of the side of a perfect cube to arrange N atoms.
    // “Ceiling” to rounding up.
    int Nx = std::ceil(std::cbrt(N));
    int Ny = Nx;
    int Nz = Nx;

    // The distribution is created here, using the external generator 'rng'.
    std::normal_distribution<double> dist(0.0, noise_scale * D);

    int count = 0;
    for (int i = 0; i < Nx && count < N; ++i) {
        for (int j = 0; j < Ny && count < N; ++j) {
            for (int k = 0; k < Nz && count < N; ++k) {
                // Use the 'rng' generator passed as an argument to create the random numbers.
                double x = i * D + dist(rng);
                double y = j * D + dist(rng);
                double z = k * D + dist(rng);
                atoms.emplace_back(x, y, z);
                ++count;
            }
        }
    }
}


std::vector<Point3D>& Cluster::getAtoms() {
    return atoms;
}

double Cluster::getEnergy() const {
    return Energy;
}




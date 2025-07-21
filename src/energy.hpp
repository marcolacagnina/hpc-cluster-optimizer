#ifndef ENERGY_HPP
#define ENERGY_HPP

#include "../include/point3d.hpp"
#include <vector>
#include <cmath>

namespace Energy {

    /**
 * @brief Calculate the total energy of a cluster of atoms using the Lennard-Jones potential.
 * This version is optimized with an OpenMP SIMD directive to vectorize
 * the computation of pair interactions, significantly speeding up execution.
 * * @param atoms 3D point vector representing the positions of atoms.
 * @return The total potential energy of the system.
 */
    inline double computeLJEnergy_SIMD(const std::vector<Point3D>& atoms) {
        const double epsilon = 1.0;
        const double sigma = 1.0;
        double energy = 0.0;
        const std::size_t N = atoms.size();

        // The constant sigma^2 is calculated only once outside of loops.
        const double sigma2 = sigma * sigma;

        // The outer loop iterates over each atom i
        for (std::size_t i = 0; i < N; ++i) {
#pragma omp simd reduction(+:energy)
            for (std::size_t j = i + 1; j < N; ++j) {
                const double dx = atoms[i].x - atoms[j].x;
                const double dy = atoms[i].y - atoms[j].y;
                const double dz = atoms[i].z - atoms[j].z;

                const double r2 = dx * dx + dy * dy + dz * dz;

                // Avoid std::pow() and std::sqrt()
                const double inv_r2 = 1.0 / r2;
                const double s2_div_r2 = sigma2 * inv_r2;
                const double r6 = s2_div_r2 * s2_div_r2 * s2_div_r2;    // (sigma^2/r^2)^3
                const double r12 = r6 * r6;

                energy += 4.0 * epsilon * (r12 - r6);
            }
        }

        return energy;
    }

} // namespace Energy

#endif

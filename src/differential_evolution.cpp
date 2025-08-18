#include "differential_evolution.hpp"
#include <random>
#include <iostream>
#include <algorithm>
#include <omp.h>


OptimizationResult differentialEvolution_Adaptive(
        int NP,
        int generations,
        int N_atoms,
        double D,
        const std::function<double(const std::vector<Point3D>&)>& energyFunc,
        double CR_fixed
) {
    // "Brutal Jump" Parameters
    const int reinitialization_stall_trigger = 100;
    const double reinitialization_fraction = 0.20;

    // --- Data Structures ---
    std::vector<Point3D> population_flat(NP * N_atoms);
    std::vector<Point3D> next_population_flat(NP * N_atoms);
    std::vector<double> energies(NP);
    std::vector<double> next_energies(NP);
    std::vector<double> F_vals(NP);

    // Centralized creations of random generators, each for each threads.
    int max_threads = omp_get_max_threads();
    std::vector<std::mt19937> rngs;
    rngs.reserve(max_threads);
    std::random_device rd;
    for (int i = 0; i < max_threads; ++i) {
        rngs.emplace_back(rd() + i); // Unique seed, each for each threads
    }

    // === Optimized Parallel Initialization ===
#pragma omp parallel
    {
        // Each thread gets its own pre-initialized generator
        auto& rng_local = rngs[omp_get_thread_num()];
        std::uniform_real_distribution<double> dist_F_range_init(0.3, 1.0);

#pragma omp for
        for (int i = 0; i < NP; ++i) {
            Cluster cl(N_atoms, D, energyFunc, rng_local);

            const auto& atoms = cl.getAtoms();
            std::copy(atoms.begin(), atoms.end(), population_flat.begin() + i * N_atoms);
            energies[i] = cl.getEnergy();
            F_vals[i] = dist_F_range_init(rng_local);
        }
    }

    // --- Serial Part before Generations ---
    double best_energy = *std::min_element(energies.begin(), energies.end());
    int stall_counter = 0;
    const int early_stopping_stall = 500;
    const double tolerance = 1e-6;

    // --- Buffer pre-allocation ---
    std::vector<std::vector<Point3D>> mutant_buffer(max_threads, std::vector<Point3D>(N_atoms));
    std::vector<std::vector<Point3D>> trial_buffer(max_threads, std::vector<Point3D>(N_atoms));

    // === Generations Principle cycle ===
    for (int gen = 0; gen < generations; ++gen) {
#pragma omp parallel
        {
            // Using local generator already pre-allocated
            int thread_id = omp_get_thread_num();
            auto& rng_local = rngs[thread_id];
            auto& mutant = mutant_buffer[thread_id];
            auto& trial = trial_buffer[thread_id];
            std::uniform_real_distribution<double> dist_real(0.0, 1.0);
            std::uniform_real_distribution<double> dist_F_range(0.3, 1.0);
#pragma omp for
            for (int i = 0; i < NP; ++i) {
                // ... INTERNAL LOGIC OF DE ...
                double F_current = F_vals[i];
                if (dist_real(rng_local) < 0.1) { F_current = dist_F_range(rng_local); }
                const double CR = CR_fixed;
                int r1, r2, r3;
                do { r1 = rng_local() % NP; } while (r1 == i);
                do { r2 = rng_local() % NP; } while (r2 == i || r2 == r1);
                do { r3 = rng_local() % NP; } while (r3 == i || r3 == r1 || r3 == r2);
                for (int j = 0; j < N_atoms; ++j) {
                    mutant[j].x = population_flat[r1 * N_atoms + j].x + F_current * (population_flat[r2 * N_atoms + j].x - population_flat[r3 * N_atoms + j].x);
                    mutant[j].y = population_flat[r1 * N_atoms + j].y + F_current * (population_flat[r2 * N_atoms + j].y - population_flat[r3 * N_atoms + j].y);
                    mutant[j].z = population_flat[r1 * N_atoms + j].z + F_current * (population_flat[r2 * N_atoms + j].z - population_flat[r3 * N_atoms + j].z);
                }
                int j_rand = rng_local() % N_atoms;
                for (int j = 0; j < N_atoms; ++j) {
                    trial[j] = (dist_real(rng_local) < CR || j == j_rand) ? mutant[j] : population_flat[i * N_atoms + j];
                }
                double trial_energy = energyFunc(trial);
                if (trial_energy < energies[i]) {
                    std::copy(trial.begin(), trial.end(), next_population_flat.begin() + i * N_atoms);
                    next_energies[i] = trial_energy;
                    F_vals[i] = F_current;
                } else {
                    std::copy(population_flat.begin() + i * N_atoms, population_flat.begin() + (i + 1) * N_atoms, next_population_flat.begin() + i * N_atoms);
                    next_energies[i] = energies[i];
                }
            }
        }

        // --- Serial part of End Generation ---
        population_flat.swap(next_population_flat);
        energies.swap(next_energies);
        double current_best_energy = *std::min_element(energies.begin(), energies.end());
        if (std::abs(current_best_energy - best_energy) < tolerance) {
            stall_counter++;
        } else {
            stall_counter = 0;
            best_energy = current_best_energy;
        }


        if (stall_counter >= early_stopping_stall) {
            if (gen < generations - 1) {
                std::cout << "Early stopping at generation " << gen
                          << " due to stagnation. Best energy: " << best_energy << std::endl;
            }
            break;          // END ALGORITHM RUN
        }


        // === RESET MECHANISM ===
        if (stall_counter > 0 && stall_counter % reinitialization_stall_trigger == 0) {
            std::cout << "  -> Generation " << gen << ": Stagnation detected. Re-initializing  "
                      << reinitialization_fraction * 100 << "% of the population." << std::endl;

            // Compute exact number of individual to reinitialize.
            const int num_to_reinitialize = static_cast<int>(NP * reinitialization_fraction);

            //Create a vector of indexes
            std::vector<int> indices(NP);
            std::iota(indices.begin(), indices.end(), 0);

            // Order indexes based on the energy of the individual they refer to.
            // At the end, vector “indices”, It will contain the indices of the best individuals (lowest energy) first
            // and those of the worst individuals (highest energy) last.
            std::stable_sort(indices.begin(), indices.end(),
                             [&](int a, int b) { return energies[a] < energies[b]; });

            // Re-initialize the worst clusters, THE BRUTAL JUMP
#pragma omp parallel for
            for (int i = 0; i < num_to_reinitialize; ++i) {

                // Select index of a cluster that is in the set of the worst cluster.
                int worst_idx = indices[NP - 1 - i];

                // Use the local random generator of the thread and pass it to the cluster constructor.
                auto& rng_local = rngs[omp_get_thread_num()];
                Cluster cl(N_atoms, D, energyFunc, rng_local);

                const auto& new_atoms = cl.getAtoms();

                // Substitute older cluster with the new already created.
                std::copy(new_atoms.begin(), new_atoms.end(), population_flat.begin() + worst_idx * N_atoms);
                energies[worst_idx] = cl.getEnergy();
            }
            best_energy = *std::min_element(energies.begin(), energies.end());
        }

    }

    // ==== FOUND FINAL RESULT AND RETURN IT ====
    auto best_it = std::min_element(energies.begin(), energies.end());
    auto best_idx = std::distance(energies.begin(), best_it);
    std::vector<Point3D> best_atoms(N_atoms);
    std::copy(population_flat.begin() + best_idx * N_atoms, population_flat.begin() + (best_idx + 1) * N_atoms, best_atoms.begin());
    return {best_atoms, *best_it};
}









#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "point3d.hpp"
#include <vector>
#include <string>
#include <functional>
#include <random>    

class Cluster {
public:

    Cluster(int N, double D, const std::function<double(const std::vector<Point3D>&)>& energyFunc, std::mt19937& rng);

    std::vector<Point3D>& getAtoms();
    double getEnergy() const;


private:
    double Energy;
    std::vector<Point3D> atoms;
    std::function<double(const std::vector<Point3D>&)> energyFunction;

    void generateRandomBoxCluster(int N, double D, std::mt19937& rng);
};

#endif // CLUSTER_HPP

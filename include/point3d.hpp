#ifndef POINT3D_HPP
#define POINT3D_HPP

#include <iostream>


class Point3D {
public:
    double x, y, z;

    Point3D();
    Point3D(double x, double y, double z);

    void print() const;
};

#endif

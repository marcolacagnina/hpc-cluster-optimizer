#include "point3d.hpp"

Point3D::Point3D() : x(0.0), y(0.0), z(0.0) {}
Point3D::Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

void Point3D::print() const {
    std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
}




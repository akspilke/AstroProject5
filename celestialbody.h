#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include "vec3.h"

class CelestialBody
{
public:
    bool bound = true;
    Vec3 position;
    Vec3 velocity;
    Vec3 force;
    double mass;

    CelestialBody(Vec3 position, Vec3 velocity, double mass);
    CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass);
    void resetForce();
    void setIsBound(bool b) { bound = b; }
    bool isBound() { return bound; }
};

#endif // CELESTIALBODY_H`

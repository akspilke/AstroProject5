#include "celestialbody.h"

CelestialBody::CelestialBody(Vec3 pos, Vec3 vel, double mass_) {
    position = pos;
    velocity = vel;
    mass = mass_;
}

CelestialBody::CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass_) {
    position = Vec3(x,y,z);
    velocity = Vec3(vx,vy,vz);
    mass = mass_;
}

void CelestialBody::resetForce() {
    force.zeros();
}

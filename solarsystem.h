#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

class SolarSystem
{
public:
    SolarSystem();
    void createCelestialBody(Vec3 position, Vec3 velocity, double mass);
    void calculateForcesAndEnergy();
    int numberOfBodies() const;
    double R0 = 20;  //light years
    double dt = .0001;
    double avg_mass = 0;
    double t_crunch = 0;

    void setTau(double tau) {
        t_crunch = tau;
    }

    void setAvgMass(double avg_m) {
        avg_mass = avg_m;
    }

    double m_G;

    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    double totalmass() const;
    double kineticBounded() const;
    double potentialBounded() const;
    double density() const;
    void calculateBoundedEnergy();
    void computeDensity();

    void writeToFile(std::string filename);
    void writeToFile2(std::string filename1, double totalK, double totalV);
    void writeToFile3(std::string , int distrubution);


    std::vector<CelestialBody> &bodies();

private:
    std::vector<CelestialBody> m_bodies;
    std::ofstream m_file;
    double m_totalmass;
    double m_kineticEnergy;
    double m_potentialEnergy;

    double m_kineticBounded   = 0;
    double m_potentialBounded = 0;
};

#endif // SOLARSYSTEM_H

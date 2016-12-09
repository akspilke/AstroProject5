#ifndef VERLET_H
#define VERLET_H


class Verlet
{
public:
    int numTimesteps = 1000;
    double m_dt;
    Verlet(double dt);
    void integrateOneStep(class SolarSystem &system);
    };


#endif // VERLET_H

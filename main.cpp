#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include "time.h"
#include <string>
#include <limits>
#include <random>
#include <armadillo>
#include "verlet.h"

using namespace std;
using namespace arma;


int main(int numArguments, char **arguments)
{
    SolarSystem solarSystem;

    ofstream outputFile;
    outputFile.open("Hist_densityN100.txt");
    outputFile << setiosflags(ios::showpoint | ios::uppercase);

    double t_crunch = 1;
    solarSystem.setTau(t_crunch);
    //double totaltime = 5*t_crunch;

    double dt = 1e-3;
    double tot_time = 5.0;
    int numTimesteps = (int) (tot_time/dt);

    if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

    //solarSystem.setTcrunch(tau);
    // We create new bodies like this. Note that the createCelestialBody function returns a reference to the newly created body
    // This can then be used to modify properties or print properties of the body if desired
    // Use with: solarSystem.createCelestialBody( position, velocity, mass );
    int N = 100;
    double R0 = 20;  //light years


    default_random_engine generator;


    //normal_distribution<double> distribution(10.0,1); //gaussian with 10Msun mean, 1Msun STDEV
    normal_distribution<double> distribution(10.0,1.0); //gaussian with 2Msun mean, 0.2Msun STDEV

    //srand(clock());
    double avg_mass = 0;
    for(int i=0; i< N;i++){
        mat randomPosition = randu(3,1);
        double u = randomPosition(0,0);
        double v = randomPosition(1,0);
        double w = randomPosition(2,0);

        //mat randomVelocity = randu(3,1);
        //double u2 = randomVelocity(0,0);
        //double v2 = randomVelocity(1,0);
        //double w2 = randomVelocity(2,0);

        double phi   = w * 2 * M_PI;
        double theta = acos(1-2*v);
        double r     = R0 * pow(u, 1./3);

        double x = r*sin(theta)*cos(phi);
        double y = r*sin(theta)*sin(phi);
        double z = r*cos(theta);

        double bodymass = distribution(generator);
        avg_mass += bodymass;

        solarSystem.createCelestialBody( Vec3(x,y,z), Vec3(0,0,0), bodymass);
    }
    avg_mass /= N;
    solarSystem.setAvgMass(avg_mass);


    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
    }

    Verlet integrator(dt);
    double totalV;
    double totalK;

    for(int timestep=0; timestep<numTimesteps; timestep++) {
        integrator.integrateOneStep(solarSystem);
        //solarSystem.writeToFile("positions_xyz_100.txt");


        if (timestep % 1 == 0) {

            solarSystem.calculateBoundedEnergy();
            solarSystem.writeToFile("Energy_unsmooth_N50.txt");

        }
    }

    solarSystem.computeDensity();
    return 0;
}


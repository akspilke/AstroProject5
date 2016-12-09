#include "solarsystem.h"
#include "vec3.h"
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <fstream>
using namespace std;




SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_totalmass(0)
{
}
void SolarSystem::createCelestialBody(Vec3 position, Vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
}
void SolarSystem::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_totalmass = 0;
    double dr = 0;

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];

        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            double dx = body1.position[0] - body2.position[0];
            double dy = body1.position[1] - body2.position[1];
            double dz = body1.position[2] - body2.position[2];
            double dr2 = dx*dx + dy*dy + dz*dz;
            double dr = sqrt(dr2);
            double epsilon = 0.1;

            m_G= (M_PI)/((32*numberOfBodies()*avg_mass*t_crunch*t_crunch)/(4*M_PI*R0*R0*R0));

            Vec3 deltaRVector = body1.position - body2.position;
            dr = deltaRVector.length();

            //code for forces on bodies without the smoothing factor
            body1.force += (-m_G*body1.mass*body2.mass*deltaRVector)/pow(dr,3);
            body2.force -= (-m_G*body1.mass*body2.mass*deltaRVector)/pow(dr,3);

            //code for forces on bodies with the smoothing factor
            //body1.force += (-m_G*body1.mass*body2.mass*deltaRVector)/(pow(dr,3)+(epsilon*epsilon));
            //body2.force -= (-m_G*body1.mass*body2.mass*deltaRVector)/(pow(dr,3)+(epsilon*epsilon));

            //m_totalmass += body1.mass+body2.mass;

            m_potentialEnergy += -(m_G*body1.mass*body2.mass)/dr;
        }
        m_totalmass     += body1.mass;
        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
}
int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}
double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}
double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}
double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

void SolarSystem::calculateBoundedEnergy() {

    int numberBounded = 0;
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &bodyi = m_bodies[i];
        bodyi.setIsBound(false);

        double energyiTotal = 0;

        // compute kinetic energy of bodyi
        energyiTotal += 0.5*bodyi.mass*bodyi.velocity.lengthSquared();

        for(int j=0; j<numberOfBodies(); j++) {
            CelestialBody &bodyj = m_bodies[j];
            if (i != j) {
                // compute potential(i,j)
                double dx = bodyi.position[0] - bodyj.position[0];
                double dy = bodyi.position[1] - bodyj.position[1];
                double dz = bodyi.position[2] - bodyj.position[2];
                double dr2 = dx*dx + dy*dy + dz*dz;
                double dr = sqrt(dr2);
                m_G= (M_PI)/((32*numberOfBodies()*avg_mass*t_crunch*t_crunch)/(4*M_PI*R0*R0*R0));
                Vec3 deltaRVector = bodyi.position - bodyj.position;
                dr = deltaRVector.length();
                energyiTotal += -(m_G*bodyi.mass*bodyj.mass)/dr; //subtract potential energy from total
            }
        }
        // Counting bounded objects
        if (energyiTotal < 0) {
            bodyi.setIsBound(true);
            numberBounded++;
        }
    }
    double totalK = 0;
    double totalV = 0;

    // Compute kinetic and potential energy only of bounded objects
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &bodyi = m_bodies[i];
        if (bodyi.isBound()) {
            // compute kinetic energy of bodyi
            totalK += 0.5*bodyi.mass*bodyi.velocity.lengthSquared();;
            for(int j=i+1; j<numberOfBodies(); j++) {
                CelestialBody &bodyj = m_bodies[j];
                // compute potential(i,j)
                if (bodyj.isBound()) {
                    double dx = bodyi.position[0] - bodyj.position[0];
                    double dy = bodyi.position[1] - bodyj.position[1];
                    double dz = bodyi.position[2] - bodyj.position[2];
                    double dr2 = dx*dx + dy*dy + dz*dz;
                    double dr = sqrt(dr2);
                    m_G= (M_PI)/((32*numberOfBodies()*avg_mass*t_crunch*t_crunch)/(4*M_PI*R0*R0*R0));
                    Vec3 deltaRVector = bodyi.position - bodyj.position;
                    dr = deltaRVector.length();
                    totalV += -(m_G*bodyi.mass*bodyj.mass)/dr;
                }
            }
        }
    }
    //writeToFile2("boundsmooth_100_e1.txt",totalK,totalV);
}

// Compute number density of particles in different shells/bins
void SolarSystem::computeDensity() {
    int numberBounded = 0;
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &bodyi = m_bodies[i];
        bodyi.setIsBound(false);
        double energyiTotal = 0;
        // compute kinetic energy of bodyi
        energyiTotal += 0.5*bodyi.mass*bodyi.velocity.lengthSquared();
        for(int j=0; j<numberOfBodies(); j++) {
            CelestialBody &bodyj = m_bodies[j];
            if (i != j) {
                // compute potential(i,j)
                double dx = bodyi.position[0] - bodyj.position[0];
                double dy = bodyi.position[1] - bodyj.position[1];
                double dz = bodyi.position[2] - bodyj.position[2];
                double dr2 = dx*dx + dy*dy + dz*dz;
                double dr = sqrt(dr2);
                m_G= (M_PI)/((32*numberOfBodies()*avg_mass*t_crunch*t_crunch)/(4*M_PI*R0*R0*R0));
                Vec3 deltaRVector = bodyi.position - bodyj.position;
                dr = deltaRVector.length();
                energyiTotal += -(m_G*bodyi.mass*bodyj.mass)/dr;
            }
        }
        if (energyiTotal < 0) {
            bodyi.setIsBound(true);
            numberBounded++;
        }
    }

    // Defining shells to go through
    double r1 = 0;          // Starting at centre of the sphere
    double r_end = 4*R0;    // Ending at 4 times the initial radius
    double step = 2;        // Defining thickness of shells, in ly
    int Nshells = (int)r_end/step;  // Total number of shells

    vector<int> distribution (Nshells);   // Initialising vector for making the histogram

    // GOING TROUGH ALL SHELLS COUNTING PARTICLES WITHIN THE SHELLS
    for(int j = 0; j < Nshells; j++) {
        r1 = j*step;
        for(int i = 0; i < numberOfBodies(); i++) {
            CelestialBody& body1 = m_bodies[i];
            if(body1.bound) {

                double dx2 = body1.position[0];
                double dy2 = body1.position[1];
                double dz2 = body1.position[2];
                double dr22 = dx2*dx2 + dy2*dy2 + dz2*dz2;
                double distanceFromCentre = sqrt(dr22);

                double ShellVol = r1*r1*4*M_PI*step;

                //cout << distanceFromCentre << endl;
                if( (distanceFromCentre >= r1 ) && ( distanceFromCentre <= r1+step) ) {
                    distribution[j]++;
                }
            }
        }
    }
    // Writing ditribution to file to make histogram
    for(int k = 0; k < Nshells; k++) {
        cout << distribution[k] << endl;
        //writeToFile3("Histogram_n2000_fixmass05.txt", distribution[k]);
    }
    //cout << numberBounded << endl;
}

//writes the positions and all energies to file
void SolarSystem::writeToFile(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }
    m_file << "100" << endl;
    m_file << "whatever" << endl;
    for(CelestialBody &body : m_bodies) {
        m_file << setprecision(10) << body.position.x() << " " << setprecision(10) << body.position.y() << " " << setprecision(10) << body.position.z() << "\n";

    }
     //m_file << setprecision(10) << totalEnergy() << " " << setprecision(10) << kineticEnergy() << " " << setprecision(10) << potentialEnergy() << "\n";
}

//write the bound energies to file
void SolarSystem::writeToFile2(string filename1, double totalK, double totalV)
{
    if(!m_file.good()) {
        m_file.open(filename1.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename1 << ". Aborting!" << endl;
            terminate();
        }
    }
   //m_file << totalV << setprecision(10) << " " << setprecision(10) << totalK << " " << "\n";
}

//write ditribution of particles to file for histogram
void SolarSystem::writeToFile3(string filename2, int distribution)
{
    if(!m_file.good()) {
        m_file.open(filename2.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename2 << ". Aborting!" << endl;
            terminate();
        }
    }
    m_file << distribution << "\n";
}

std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}

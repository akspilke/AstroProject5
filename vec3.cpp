#include "vec3.h"
#include <cmath>
#include <iostream>
#include <armadillo>
using namespace std;


Vec3::Vec3()
{
    zeros();
}

Vec3::Vec3(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

void Vec3::print()
{
    // Will print matlab syntax vector. Output will be like: [2.09, 5.3, 9.1];
    cout << "[" << components[0] << ", " << components[1] << ", " << components[2] << "]" << endl;
}

void Vec3::print(string name)
{
    // Will print matlab syntax vector with a name. Output will be like: A = [2.09, 5.3, 9.1];
    cout << name << " = ";
    print();
}

Vec3 Vec3::cross(Vec3 otherVector)
{
    return Vec3(y()*otherVector.z()-z()*otherVector.y(), z()*otherVector.x()-x()*otherVector.z(), x()*otherVector.y()-y()*otherVector.x());
}

double Vec3::dot(Vec3 otherVector)
{
    return otherVector[0]*components[0] + otherVector[1]*components[1] + otherVector[2]*components[2];
}

void Vec3::normalize()
{
    double length = this->length();
    if(length > 0) {
        components[0] /= length;
        components[1] /= length;
        components[2] /= length;
    }
}

Vec3 Vec3::normalized()
{
    Vec3 newVector = *this;
    newVector.normalize();
    return newVector;
}

double Vec3::lengthSquared()
{
    // Returns the square of the length (or norm) of the vector
    return components[0]*components[0]+components[1]*components[1]+components[2]*components[2];
}

double Vec3::length()
{
    // Returns the length (or norm) of the vector
    return sqrt(lengthSquared());
}

void Vec3::zeros()
{
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
}

Vec3 &Vec3::operator+=(double rhs)
{
    components[0] += rhs;
    components[1] += rhs;
    components[2] += rhs;
    return *this;
}

Vec3 &Vec3::operator+=(Vec3 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    components[2] += rhs[2];
    return *this;
}

Vec3 &Vec3::operator*=(double rhs)
{
    components[0] *= rhs;
    components[1] *= rhs;
    components[2] *= rhs;
    return *this;
}

Vec3 &Vec3::operator*=(Vec3 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    components[2] *= rhs[2];
    return *this;
}

Vec3 &Vec3::operator-=(double rhs)
{
    components[0] -= rhs;
    components[1] -= rhs;
    components[2] -= rhs;
    return *this;
}

Vec3 &Vec3::operator-=(Vec3 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    components[2] -= rhs[2];
    return *this;
}

Vec3 &Vec3::operator/=(double rhs)
{
    components[0] /= rhs;
    components[1] /= rhs;
    components[2] /= rhs;
    return *this;
}

Vec3 &Vec3::operator/=(Vec3 rhs)
{
    components[0] /= rhs[0];
    components[1] /= rhs[1];
    components[2] /= rhs[2];
    return *this;
}

std::ostream &operator<<(std::ostream &os, const Vec3 &myVector) // Allows cout << myVector << endl;
{
    os << "[" << myVector.x() << ", " << myVector.y() << ", " << myVector.z() << "];";
    return os;
}

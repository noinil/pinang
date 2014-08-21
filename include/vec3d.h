// -*-c++-*-

#ifndef PINANG_VEC3D_H_
#define PINANG_VEC3D_H_

#include "constants.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

namespace pinang{
    class Vec3d
    {
    public:
        Vec3d(): _z1(0), _z2(0), _z3(0) {}
        Vec3d(double a, double b, double c): _z1(a), _z2(b), _z3(c) {}

        double x() const {return _z1;}
        double y() const {return _z2;}
        double z() const {return _z3;}

        Vec3d operator+(const Vec3d& other) const
        {
            return Vec3d(_z1+other._z1,_z2+other._z2,_z3+other._z3);
        }

        Vec3d operator-(const Vec3d& other) const
        {
            return Vec3d(_z1-other._z1,_z2-other._z2,_z3-other._z3);
        }

        Vec3d operator*(const double x) const
        {
            return Vec3d(_z1*x,_z2*x,_z3*x);
        }

        double operator*(const Vec3d& other) const
        {
            return (_z1*other._z1+_z2*other._z2+_z3*other._z3);
        }

        Vec3d operator^(const Vec3d& other) const
        {
            return Vec3d(_z2*other._z3 - _z3*other._z2,
                         _z3*other._z1 - _z1*other._z3,
                         _z1*other._z2 - _z2*other._z1);
        }

        double operator[](unsigned int i) const
        {
            switch(i)
            {
            case 0: return _z1;
            case 1: return _z2;
            case 2: return _z3;
            default:
                std::cerr << "ERROR: 3D Vector index out of range!" << std::endl;
                exit(EXIT_SUCCESS);
            }
        }

        inline double norm() const
        {
            return sqrt(_z1*_z1+_z2*_z2+_z3*_z3);
        }

        inline double sqr_norm() const
        {
            return (_z1*_z1+_z2*_z2+_z3*_z3);
        }

    protected:
        double _z1, _z2, _z3;
    };

    inline std::ostream& operator<<(std::ostream& o, const Vec3d& v)
    {
        o << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
          << std::setw(8) << v.x()
          << std::setw(8) << v.y()
          << std::setw(8) << v.z();
        // std::cout.unsetf(std::ios::floatfield);
        // std::cout.precision(6);
        return o;
    }

    inline std::istream& operator>>(std::istream& i, Vec3d& v)
    {
        double x,y,z;
        i >> x >> y >> z;
        if (!i) return i;
        v = Vec3d(x,y,z);
        return i;
    }

    inline double vec_distance (Vec3d& v1, Vec3d& v2)
    {
        // Actually this is not the distance between two Vectors,
        // but the distance between two Points!
        Vec3d v3 = v1-v2;
        return v3.norm();
    }

    double vec_angle (Vec3d& v1, Vec3d& v2)
    {
        if ((v1*v2)/(v1.norm()*v2.norm()) >= 0.99999)
            return 0;
        if ((v1*v2)/(v1.norm()*v2.norm()) <= -0.99999)
            return 3.1415926;
        return acos((v1*v2)/(v1.norm()*v2.norm()));
    }
    double vec_angle_deg (Vec3d& v1, Vec3d& v2)
    {
        if ((v1*v2)/(v1.norm()*v2.norm()) >= 0.99999)
            return 0;
        if ((v1*v2)/(v1.norm()*v2.norm()) <= -0.99999)
            return 180;
        double ang = acos((v1*v2)/(v1.norm()*v2.norm()));
        return (180 * ang / g_pi);
    }

}
#endif

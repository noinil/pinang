// -*-c++-*-

#ifndef BINANG_VEC3D_H
#define BINANG_VEC3D_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

namespace binang{
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
                assert(i < 3);
                return -1;
            }
        }

        double norm() const
        {
            return sqrt(_z1*_z1+_z2*_z2+_z3*_z3);
        }
    protected:
        double _z1, _z2, _z3;
    };

    inline std::ostream& operator<<(std::ostream& o, const Vec3d& v)
    {
        o << std::setw(5) << v.x() << " "
          << std::setw(5) << v.y() << " "
          << std::setw(5) << v.z();
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

    inline double Vec_distance (Vec3d& v1, Vec3d& v2)
    {
        Vec3d v3 = v1-v2;
        return v3.norm();
    }

}
#endif

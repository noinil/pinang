// -*-c++-*-

#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

namespace binang{
    class Vec3d
    {
    private:
        double z1, z2, z3;
    public:
        Vec3d(): z1(0), z2(0), z3(0) {}
        Vec3d(double a, double b, double c): z1(a), z2(b), z3(c) {}

        double x() const {return z1;}
        double y() const {return z2;}
        double z() const {return z3;}

        Vec3d operator+(const Vec3d &other) const
        {
            return Vec3d(z1+other.z1,z2+other.z2,z3+other.z3);
        }

        Vec3d operator-(const Vec3d &other) const
        {
            return Vec3d(z1-other.z1,z2-other.z2,z3-other.z3);
        }

        Vec3d operator*(const double x) const
        {
            return Vec3d(z1*x,z2*x,z3*x);
        }

        double operator*(const Vec3d &other) const
        {
            return (z1*other.z1+z2*other.z2+z3*other.z3);
        }

        Vec3d operator^(const Vec3d &other) const
        {
            return Vec3d(z2*other.z3-z3*other.z2, z3*other.z1 - z1*other.z3, z1*other.z2 - z2*other.z1);
        }

        double operator[](unsigned int i) const
        {
            switch(i)
            {
            case 0: return z1;
            case 1: return z2;
            case 2: return z3;
            default:
                assert(i < 3);
                return -1;
            }
        }

        double norm() const
        {
            return sqrt(z1*z1+z2*z2+z3*z3);
        }
    };

    inline std::ostream &operator<<(std::ostream &o, const Vec3d &v)
    {
        o << std::setw(5) << v.x() << "  "
          << std::setw(5) << v.y() << "  "
          << std::setw(5) << v.z();
        return o;
    }

    inline std::istream &operator>>(std::istream &i, Vec3d &v)
    {
        double x,y,z;
        i >> x >> y >> z;
        if (!i) return i;
        v = Vec3d(x,y,z);
        return i;
    }

    inline double Vec_distance (Vec3d &v1, Vec3d &v2)
    {
        Vec3d v3 = v1-v2;
        return v3.norm();
    }

}

/*!
@file vec3d.cpp
@brief Define functions of class Vec3d.

Definitions of member or friend functions of class Vec3d.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-24 15:47
@copyright GNU Public License V3.0
*/


#include <cmath>
#include <iomanip>
#include "vec3d.hpp"

namespace pinang {

// ------------------------------ MEMBERS --------------------
double Vec3d::operator[](unsigned int i) const
{
  switch(i)
  {
    case 0: return z1_;
    case 1: return z2_;
    case 2: return z3_;
    default:
      std::cout << " ~             PINANG :: vec3d.hpp              ~ " << "\n";
      std::cerr << "ERROR: 3D Vector index out of range!" << "\n";
      exit(EXIT_SUCCESS);
  }
}

Vec3d Vec3d::operator-() const
{
  Vec3d v0;
  v0.z1_ = - z1_;
  v0.z2_ = - z2_;
  v0.z3_ = - z3_;
  return v0;
}

Vec3d Vec3d::unitv() const
{
  Vec3d v0;
  double d = sqrt(z1_ * z1_ + z2_ * z2_ + z3_ * z3_);
  if (d < 0.00000001) {
    std::cout << " ~             PINANG :: vec3d.hpp              ~ " << "\n";
    std::cerr << "ERROR: Vector length ~ 0!" << "\n";
    exit(EXIT_SUCCESS);
  }
  v0.z1_ = z1_ / d;
  v0.z2_ = z2_ / d;
  v0.z3_ = z3_ / d;
  return v0;
}
double Vec3d::norm() const
{
  double d = sqrt(z1_ * z1_ + z2_ * z2_ + z3_ * z3_);
  return d;
}

double Vec3d::squared_norm() const
{
  double d = z1_ * z1_ + z2_ * z2_ + z3_ * z3_;
  return d;
}

// ------------------------------ FRIENDS --------------------
Vec3d operator+(const Vec3d& v1, const Vec3d& v2)
{
  return Vec3d(v1.z1_ + v2.z1_, v1.z2_ + v2.z2_, v1.z3_ + v2.z3_);
}

Vec3d operator-(const Vec3d& v1, const Vec3d& v2)
{
  return Vec3d(v1.z1_ - v2.z1_, v1.z2_ - v2.z2_, v1.z3_ - v2.z3_);
}

Vec3d operator*(const Vec3d & v, double x)
{
  return Vec3d(v.z1_ * x, v.z2_ * x, v.z3_ * x);
}

Vec3d operator/(const Vec3d & v, double x)
{
  return Vec3d(v.z1_ / x, v.z2_ / x, v.z3_ / x);
}

Vec3d operator*(double x, const Vec3d & v)
{
  return Vec3d(v.z1_ * x, v.z2_ * x, v.z3_ * x);
}

double operator*(const Vec3d& v1, const Vec3d& v2)
{
  double d = v1.z1_ * v2.z1_ + v1.z2_ * v2.z2_ + v1.z3_ * v2.z3_;
  return d;
}

Vec3d operator%(const Vec3d& u, const Vec3d& v)
{
  return Vec3d(u.z2_ * v.z3_ - u.z3_ * v.z2_,
               u.z3_ * v.z1_ - u.z1_ * v.z3_,
               u.z1_ * v.z2_ - u.z2_ * v.z1_);
}

Vec3d& operator+=(Vec3d& v0, const Vec3d& v1)
{
  v0.z1_ += v1.z1_;
  v0.z2_ += v1.z2_;
  v0.z3_ += v1.z3_;
  return v0;
}
Vec3d& operator-=(Vec3d& v0, const Vec3d& v1)
{
  v0.z1_ -= v1.z1_;
  v0.z2_ -= v1.z2_;
  v0.z3_ -= v1.z3_;
  return v0;
}
Vec3d& operator*=(Vec3d& v0, double d1)
{
  v0.z1_ *= d1;
  v0.z2_ *= d1;
  v0.z3_ *= d1;
  return v0;
}
Vec3d& operator/=(Vec3d& v0, double d1)
{
  v0.z1_ /= d1;
  v0.z2_ /= d1;
  v0.z3_ /= d1;
  return v0;
}
// ---------- other --------------------
std::ostream& operator<<(std::ostream& o, const Vec3d& v)
{
  o << std::setiosflags(std::ios_base::fixed)
    << std::setprecision(6)
    << std::setw(16) << v.z1_ << " "
    << std::setw(16) << v.z2_ << " "
    << std::setw(16) << v.z3_;
  return o;
}

std::istream& operator>>(std::istream& i, Vec3d& v)
{
  double x, y, z;
  i >> x >> y >> z;
  if (!i) return i;
  v.z1_ = x;
  v.z2_ = y;
  v.z3_ = z;
  return i;
}

double vec_distance(const Vec3d& v1, const Vec3d& v2)
{
  Vec3d v3 = v1-v2;
  double d = v3.norm();
  return d;
}

double vec_angle(const Vec3d& v1, const Vec3d& v2)
{
  double d = v1 * v2;
  double m = v1.norm() * v2.norm();
  double f = d / m;
  if (f >= 0.99999)
    return 0.0;
  if (f <= -0.99999)
    return 3.1415926;
  double c = acos(f);
  return c;
}

double vec_angle_deg(const Vec3d& v1, const Vec3d& v2)
{
  double d = v1 * v2;
  double m = v1.norm() * v2.norm();
  double f = d / m;
  if (f >= 0.99999)
    return 0.0;
  if (f <= -0.99999)
    return 180.0;
  double ang = acos(f);
  d = (180.0 * ang / 3.1415926535897933);
  return d;
}

}  // pinang

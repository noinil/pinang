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
  Vec3d(): z1_(0), z2_(0), z3_(0) {}
  Vec3d(double a, double b, double c): z1_(a), z2_(b), z3_(c) {}

  inline double x() const {return z1_;}
  inline double y() const {return z2_;}
  inline double z() const {return z3_;}

  inline Vec3d operator+(const Vec3d&) const;
  inline Vec3d operator-(const Vec3d&) const;
  inline Vec3d operator*(double) const;
  inline double operator*(const Vec3d&) const;
  inline Vec3d operator^(const Vec3d&) const;
  inline double operator[](unsigned int) const;

  inline double norm() const;
  inline double sqr_norm() const;

  friend inline Vec3d operator*(double, const Vec3d&);
  friend inline std::ostream& operator<<(std::ostream&, const Vec3d&);
  friend inline std::istream& operator>>(std::istream&, Vec3d&);
  friend inline double vec_distance(Vec3d&, Vec3d&);
  friend inline double vec_angle(Vec3d&, Vec3d&);
  friend inline double vec_angle_deg(Vec3d&, Vec3d&);

 protected:
  double z1_, z2_, z3_;
};
// ------------------------------ MEMBERS --------------------
inline Vec3d Vec3d::operator+(const Vec3d& other) const
{
  return Vec3d(z1_ + other.z1_, z2_ + other.z2_, z3_ + other.z3_);
}

inline Vec3d Vec3d::operator-(const Vec3d& other) const
{
  return Vec3d(z1_ - other.z1_, z2_ - other.z2_, z3_ - other.z3_);
}

inline Vec3d Vec3d::operator*(double x) const
{
  return Vec3d(z1_ * x, z2_ * x, z3_ * x);
}

inline double Vec3d::operator*(const Vec3d& other) const
{
  return (z1_ * other.z1_ + z2_ * other.z2_ + z3_ * other.z3_);
}

inline Vec3d Vec3d::operator^(const Vec3d& other) const
{
  return Vec3d(z2_ * other.z3_ - z3_ * other.z2_,
               z3_ * other.z1_ - z1_ * other.z3_,
               z1_ * other.z2_ - z2_ * other.z1_);
}

inline double Vec3d::operator[](unsigned int i) const
{
  switch(i)
  {
    case 0: return z1_;
    case 1: return z2_;
    case 2: return z3_;
    default:
      std::cerr << "ERROR: 3D Vector index out of range!" << std::endl;
      exit(EXIT_SUCCESS);
  }
}

inline double Vec3d::norm() const
{
  return sqrt(z1_ * z1_ + z2_ * z2_ + z3_ * z3_);
}

inline double Vec3d::sqr_norm() const
{
  return (z1_ * z1_ + z2_ * z2_ + z3_ * z3_);
}

// ------------------------------ FRIENDS --------------------
inline Vec3d operator*(double x, const Vec3d & v)
{
  return Vec3d(v.z1_ * x, v.z2_ * x, v.z3_ * x);
}

inline std::ostream& operator<<(std::ostream& o, const Vec3d& v)
{
  o << std::setiosflags(std::ios_base::fixed)
    << std::setprecision(6)
    << std::setw(16) << v.x() << " "
    << std::setw(16) << v.y() << " "
    << std::setw(16) << v.z();
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

inline double vec_distance(Vec3d& v1, Vec3d& v2)
{
  Vec3d v3 = v1-v2;
  return v3.norm();
}

inline double vec_angle(Vec3d& v1, Vec3d& v2)
{
  double d = v1 * v2;
  double m = v1.norm() * v2.norm();
  double f = d / m;
  if (f >= 0.99999)
    return 0.0;
  if (f <= -0.99999)
    return k_pi;
  return acos(f);
}

inline double vec_angle_deg(Vec3d& v1, Vec3d& v2)
{
  double d = v1 * v2;
  double m = v1.norm() * v2.norm();
  double f = d / m;
  if (f >= 0.99999)
    return 0.0;
  if (f <= -0.99999)
    return 180.0;
  double ang = acos(f);
  return (180 * ang / k_pi);
}

}
#endif

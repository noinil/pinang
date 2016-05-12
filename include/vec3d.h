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

  inline double operator[](unsigned int) const;

  inline double norm() const;
  inline double sqr_norm() const;

  friend inline Vec3d operator+(const Vec3d&, const Vec3d&);
  friend inline Vec3d operator-(const Vec3d&, const Vec3d&);
  friend inline Vec3d operator*(const Vec3d&, double);  // scalar multiplication;
  friend inline Vec3d operator*(double, const Vec3d&);  // scalar multiplication;
  friend inline Vec3d operator/(const Vec3d&, double);
  friend inline Vec3d operator%(const Vec3d&, const Vec3d&);  // cross product;
  friend inline double operator*(const Vec3d&, const Vec3d&);  // dot product;

  friend inline Vec3d& operator+=(Vec3d&, const Vec3d&);
  friend inline Vec3d& operator-=(Vec3d&, const Vec3d&);
  friend inline Vec3d& operator*=(Vec3d&, double);
  friend inline Vec3d& operator/=(Vec3d&, double);

  friend inline std::ostream& operator<<(std::ostream&, const Vec3d&);
  friend inline std::istream& operator>>(std::istream&, Vec3d&);

  friend inline double vec_distance(const Vec3d&, const Vec3d&);
  friend inline double vec_angle(const Vec3d&, const Vec3d&);
  friend inline double vec_angle_deg(const Vec3d&, const Vec3d&);

 protected:
  double z1_, z2_, z3_;
};
// ------------------------------ MEMBERS --------------------
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
  double d = sqrt(z1_ * z1_ + z2_ * z2_ + z3_ * z3_);
  return d;
}

inline double Vec3d::sqr_norm() const
{
  double d = z1_ * z1_ + z2_ * z2_ + z3_ * z3_;
  return d;
}

// ------------------------------ FRIENDS --------------------
inline Vec3d operator+(const Vec3d& v1, const Vec3d& v2)
{
  return Vec3d(v1.z1_ + v2.z1_, v1.z2_ + v2.z2_, v1.z3_ + v2.z3_);
}

inline Vec3d operator-(const Vec3d& v1, const Vec3d& v2)
{
  return Vec3d(v1.z1_ - v2.z1_, v1.z2_ - v2.z2_, v1.z3_ - v2.z3_);
}

inline Vec3d operator*(const Vec3d & v, double x)
{
  return Vec3d(v.z1_ * x, v.z2_ * x, v.z3_ * x);
}

inline Vec3d operator/(const Vec3d & v, double x)
{
  return Vec3d(v.z1_ / x, v.z2_ / x, v.z3_ / x);
}

inline Vec3d operator*(double x, const Vec3d & v)
{
  return Vec3d(v.z1_ * x, v.z2_ * x, v.z3_ * x);
}

inline double operator*(const Vec3d& v1, const Vec3d& v2)
{
  double d = v1.z1_ * v2.z1_ + v1.z2_ * v2.z2_ + v1.z3_ * v2.z3_;
  return d;
}

inline Vec3d operator%(const Vec3d& u, const Vec3d& v)
{
  return Vec3d(u.z2_ * v.z3_ - u.z3_ * v.z2_,
               u.z3_ * v.z1_ - u.z1_ * v.z3_,
               u.z1_ * v.z2_ - u.z2_ * v.z1_);
}

inline Vec3d& operator+=(Vec3d& v0, const Vec3d& v1)
{
  v0.z1_ += v1.z1_;
  v0.z2_ += v1.z2_;
  v0.z3_ += v1.z3_;
  return v0;
}
inline Vec3d& operator-=(Vec3d& v0, const Vec3d& v1)
{
  v0.z1_ -= v1.z1_;
  v0.z2_ -= v1.z2_;
  v0.z3_ -= v1.z3_;
  return v0;
}
inline Vec3d& operator*=(Vec3d& v0, double d1)
{
  v0.z1_ *= d1;
  v0.z2_ *= d1;
  v0.z3_ *= d1;
  return v0;
}
inline Vec3d& operator/=(Vec3d& v0, double d1)
{
  v0.z1_ /= d1;
  v0.z2_ /= d1;
  v0.z3_ /= d1;
  return v0;
}
// ---------- other --------------------
inline std::ostream& operator<<(std::ostream& o, const Vec3d& v)
{
  o << std::setiosflags(std::ios_base::fixed)
    << std::setprecision(6)
    << std::setw(16) << v.z1_ << " , "
    << std::setw(16) << v.z2_ << " , "
    << std::setw(16) << v.z3_;
  return o;
}

inline std::istream& operator>>(std::istream& i, Vec3d& v)
{
  double x, y, z;
  i >> x >> y >> z;
  if (!i) return i;
  v.z1_ = x;
  v.z2_ = y;
  v.z3_ = z;
  return i;
}

inline double vec_distance(const Vec3d& v1, const Vec3d& v2)
{
  Vec3d v3 = v1-v2;
  double d = v3.norm();
  return d;
}

inline double vec_angle(const Vec3d& v1, const Vec3d& v2)
{
  double d = v1 * v2;
  double m = v1.norm() * v2.norm();
  double f = d / m;
  if (f >= 0.99999)
    return 0.0;
  if (f <= -0.99999)
    return k_pi;
  double c = acos(f);
  return c;
}

inline double vec_angle_deg(const Vec3d& v1, const Vec3d& v2)
{
  double d = v1 * v2;
  double m = v1.norm() * v2.norm();
  double f = d / m;
  if (f >= 0.99999)
    return 0.0;
  if (f <= -0.99999)
    return 180.0;
  double ang = acos(f);
  d = (180 * ang / k_pi);
  return d;
}

}
#endif

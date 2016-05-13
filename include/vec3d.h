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

  double x() const {return z1_;}
  double y() const {return z2_;}
  double z() const {return z3_;}

  double operator[](unsigned int) const;

  void set_coords(double x, double y, double z)
  {
    z1_ = x;
    z2_ = y;
    z3_ = z;
  }
  double norm() const;
  double sqr_norm() const;
  Vec3d unitv() const;

  Vec3d operator-() const;

  friend Vec3d operator+(const Vec3d&, const Vec3d&);
  friend Vec3d operator-(const Vec3d&, const Vec3d&);
  friend Vec3d operator*(const Vec3d&, double);  // scalar multiplication;
  friend Vec3d operator*(double, const Vec3d&);  // scalar multiplication;
  friend Vec3d operator/(const Vec3d&, double);
  friend Vec3d operator%(const Vec3d&, const Vec3d&);  // cross product;
  friend double operator*(const Vec3d&, const Vec3d&);  // dot product;

  friend Vec3d& operator+=(Vec3d&, const Vec3d&);
  friend Vec3d& operator-=(Vec3d&, const Vec3d&);
  friend Vec3d& operator*=(Vec3d&, double);
  friend Vec3d& operator/=(Vec3d&, double);

  friend std::ostream& operator<<(std::ostream&, const Vec3d&);
  friend std::istream& operator>>(std::istream&, Vec3d&);

  friend double vec_distance(const Vec3d&, const Vec3d&);
  friend double vec_angle(const Vec3d&, const Vec3d&);
  friend double vec_angle_deg(const Vec3d&, const Vec3d&);

 protected:
  double z1_, z2_, z3_;
};

}

#endif

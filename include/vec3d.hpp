/*!
@file vec3d.hpp
@brief Definition of class Vec3d.

In this file the most basic class in PINANG, Vec3d, is defined.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:12
@copyright GNU Public License V3.0
*/

#ifndef PINANG_VEC3D_H_
#define PINANG_VEC3D_H_

#include <iostream>

namespace pinang{

/*!
@brief The most basic concept of 3D coordinates.

Consisting of three components (x, y, z), the class Vec3d can perform many types
of arithmetic computations, including cross/dot products.
*/
class Vec3d
{
 public:
  //! @brief Create a Vector, whose three components all equal to 0.0.
  //! @return A Vec3d vector.
  Vec3d(): z1_(0.0), z2_(0.0), z3_(0.0) {}
  //! @brief Create a Vector based on three real numbers.
  //! @return A Vec3d vector.
  Vec3d(double a, double b, double c): z1_(a), z2_(b), z3_(c) {}

  //! @brief Get the x coordinate of a Vector.
  double x() const {return z1_;}
  //! @brief Get the y coordinate of a Vector.
  double y() const {return z2_;}
  //! @brief Get the z coordinate of a Vector.
  double z() const {return z3_;}

  //! @brief Get a component of a Vector.
  //! @param 0: x.
  //! @param 1: y.
  //! @param 2: z.
  double operator[](unsigned int) const;

  //! @brief Set the value of a Vector based on three real numbers.
  void set_coordinate(double x, double y, double z)
  {
    z1_ = x;
    z2_ = y;
    z3_ = z;
  }
  //! @brief Get the Norm of a vector.
  //! @return The Norm of a vector.
  double norm() const;
  //! @brief Get the Squared Norm of a vector.
  //! @return The Squared Norm of a vector.
  double squared_norm() const;
  Vec3d unitv() const;

  Vec3d operator-() const;

  //! @brief Perform the arithmetic computation + for two vectors.
  friend Vec3d operator+(const Vec3d&, const Vec3d&);
  //! @brief Perform the arithmetic computation - for two vectors.
  friend Vec3d operator-(const Vec3d&, const Vec3d&);
  //! @brief Perform scalar multiplication for a vector.
  friend Vec3d operator*(const Vec3d&, double);  // scalar multiplication;
  //! @brief Perform scalar multiplication for a vector.
  friend Vec3d operator*(double, const Vec3d&);  // scalar multiplication;
  //! @brief Perform scalar division for a vector.
  friend Vec3d operator/(const Vec3d&, double);
  //! @brief Get the cross product for two vectors.
  //! @return The cross product for two vectors.
  friend Vec3d operator%(const Vec3d&, const Vec3d&);  // cross product;
  //! @brief Get the dot product for two vectors.
  //! @return The dot product for two vectors.
  friend double operator*(const Vec3d&, const Vec3d&);  // dot product;

  //! @brief Add another vector to a vector.
  friend Vec3d& operator+=(Vec3d&, const Vec3d&);
  //! @brief Substract another vector from a vector.
  friend Vec3d& operator-=(Vec3d&, const Vec3d&);
  //! @brief Perform scalar multiplication for a vector.
  friend Vec3d& operator*=(Vec3d&, double);
  //! @brief Perform scalar division for a vector.
  friend Vec3d& operator/=(Vec3d&, double);

  //! @brief Output Vec3d vector.
  friend std::ostream& operator<<(std::ostream&, const Vec3d&);
  //! @brief Input Vec3d vector.
  friend std::istream& operator>>(std::istream&, Vec3d&);

  //! @brief Compute distance between two Vectors (points).
  //! @param Two Vec3d vectors.
  //! @return Distance.
  friend double vec_distance(const Vec3d&, const Vec3d&);
  //! @brief Compute angle formed between two Vectors.
  //! @param Two Vec3d vectors.
  //! @return Angle.
  friend double vec_angle(const Vec3d&, const Vec3d&);
  //! @brief Compute angle formed between two Vectors (in the unit of degree).
  //! @param Two Vec3d vectors.
  //! @return Angle (in the unit of degree).
  friend double vec_angle_deg(const Vec3d&, const Vec3d&);

 protected:
  double z1_;
  double z2_;
  double z3_;
};

double vec_distance(const Vec3d&, const Vec3d&);
double vec_angle(const Vec3d&, const Vec3d&);
double vec_angle_deg(const Vec3d&, const Vec3d&);
}

#endif

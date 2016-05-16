/*!
************************************************************
@file geometry.hpp
@brief Definition of geometric concepts.

Class Quaternion is a type to store the (mathematical) quaternion which can be
used to represent rotation matrix.  Class Transform is a class of geometrical
transforms (including rotation and translation).

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:00
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_GEOMETRY_H
#define PINANG_GEOMETRY_H

#include "vec3d.hpp"
#include "group.hpp"

namespace pinang {

class Group;

class Quaternion
{
 public:
  Quaternion(): w_(0.0), x_(0.0), y_(0.0), z_(0.0) {}
  Quaternion(double w, double x, double y, double z): w_(w), x_(x), y_(y), z_(z) {}

  double w() const {return w_;}
  double x() const {return x_;}
  double y() const {return y_;}
  double z() const {return z_;}
  int normalize();

 protected:
  double w_;  //!< The first element of a quaternion (w).
  double x_;  //!< The first component of the vector part of a quaternion (x).
  double y_;  //!< The second component of the vector part of a quaternion (y).
  double z_;  //!< The third component of the vector part of a quaternion (z).
};

class Transform
{
 public:
  Transform();
  Transform(const Quaternion&, const Vec3d&);

  Quaternion rotation() { return rotation_; }
  Vec3d translation() { return translation_; }
  void set_rotation(const Quaternion&);
  void set_translation(const Vec3d& v) { translation_ = v; }

  Vec3d apply(const Vec3d&);
  Group apply(const Group&);

 protected:
  Quaternion rotation_;  //!< The "linear" part (rotation) of a transform.
  Vec3d translation_;    //!< The translation part of a transform.

  // three vec3d for rotation matrix;
  Vec3d rotv1_;  //!< The first row of rotation matrix converted from quaternion.
  Vec3d rotv2_;  //!< The second row of rotation matrix converted from quaternion.
  Vec3d rotv3_;  //!< The third row of rotation matrix converted from quaternion.
};

}

#endif

/*!
@file geometry.hpp
@brief Definition of geometric concepts.

Class Quaternion is a type to store the (mathematical) quaternion which can be
used to represent rotation matrix.  Class Transform is a class of geometrical
transforms (including rotation and translation).

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:00
@copyright GNU Public License V3.0
*/

#ifndef PINANG_GEOMETRY_H
#define PINANG_GEOMETRY_H

#include "vec3d.hpp"
#include "group.hpp"

namespace pinang {

class Group;

/*!
@brief Simply a set of four real numbers.

In this stage, class quaternion does not follow its mathematical definition.
But I used it to denote rotation operations in geometry, mainly used in class
Transform.
*/
class Quaternion
{
 public:
  //! @brief Create a Quaternion object with all components equals to 0.
  Quaternion(): w_(0.0), x_(0.0), y_(0.0), z_(0.0) {}
  //! @brief Create a Quaternion object from four real numbers.
  Quaternion(double w, double x, double y, double z): w_(w), x_(x), y_(y), z_(z) {}

  //! @brief Get the first element of a Quaternion.
  double w() const {return w_;}
  //! @brief Get the first component of the vector part of a quaternion (x).
  double x() const {return x_;}
  //! @brief Get the second component of the vector part of a quaternion (x).
  double y() const {return y_;}
  //! @brief Get the third component of the vector part of a quaternion (x).
  double z() const {return z_;}
  //! @brief Normalize the Quaternion.
  //! @return Status of normalization.
  int normalize();

 protected:
  double w_;  //!< The first element of a quaternion (w).
  double x_;  //!< The first component of the vector part of a quaternion (x).
  double y_;  //!< The second component of the vector part of a quaternion (y).
  double z_;  //!< The third component of the vector part of a quaternion (z).
};

/*!
@brief Geometric affine operation.

Including two parts:
- The linear part: rotation.  Rrepresented by a quaternion.
- The translation part: represented by a Vec3d vector.
*/
class Transform
{
 public:
  //! @brief Create an "empty" Transform object.
  //! @return An Transform object.
  Transform();
  //! @brief Create a Transform object from a quaternion and a Vec3d type vector.
  //! @return An Transform object.
  Transform(const Quaternion&, const Vec3d&);

  //! @brief Get the rotation part of the Quaternion.
  //! @return The rotation part of the Quaternion.
  Quaternion rotation() { return rotation_; }
  //! @brief Get the translation part of the Quaternion.
  //! @return The translation part of the Quaternion.
  Vec3d translation() { return translation_; }
  //! @brief Set the rotation part of Transform object base on a quaternion.
  void set_rotation(const Quaternion&);
  //! @brief Set the translation part of Transform object base on a vector.
  void set_translation(const Vec3d& v) { translation_ = v; }

  //! @brief Apply the transform to a vector. (Create a new one)
  //! @return Vec3d which is a transformation of the original one.
  Vec3d apply(const Vec3d&);
  //! @brief Apply the transform to a Group. (Create a new one)
  //! @return Conformation which is a transformation of the original one.
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

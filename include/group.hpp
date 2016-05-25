/*!
  @file group.hpp
  @brief Definition of class group.

  In this file class Group is defined, which is basically same as class
  Conformation, but contains more member and friend functions.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 18:03
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_GROUP_H
#define PINANG_GROUP_H

#include "geometry.hpp"
#include "conformation.hpp"

namespace pinang {


class Transform;

/*!
  @brief A collection or sub-part of a conformation.

  A group of coordinates, usually as a sub-group of a conformation.  Geometric or
  Arithmetic claculations can be computed for different groups, such as RMSD,
  superimposition, or center-of-mass distances.
*/
class Group : public Conformation
{
 public:
  //! @brief Create an "empty" Group object.
  //! @return A Group object.
  Group();
  //! @brief Create a Group object based on a set of coordinates.
  //! @param A set of Vec3d coordinates.
  //! @return A Group object.
  Group(std::vector<Vec3d>);
  virtual ~Group() {coordinates_.clear();}
  //! @brief Set Group conformation based on a set of coordinates.
  //! @param A set of Vec3d coordinates, whose size should be same as original Group.
  //! @return Status of setting up the Group.
  //! @retval 1: Failure.
  //! @retval 0: Success.
  int set_conformation(std::vector<Vec3d>);

  //! @brief Get the (geometric) centroid of the Group.
  //! @return Vec3d type coordinate centroid of the Group.
  Vec3d get_centroid() const;

  //! @brief Get the center-of-mass of Group.
  //! @param Group
  //! @param A list of masses.
  //! @return Vec3d type coordinate centroid of the COM.
  friend Vec3d get_center_of_mass(const Group&, std::vector<double>);
  //! @brief Get the radius of gyration of a Group.
  //! @param Group
  //! @return Radius of gyration of a Group.
  friend double get_radius_of_gyration(const Group&);
  //! @brief Get the transform matrix from one Group to another.
  //! @param Two Group's.
  //! @param A Transform object to be calculated.
  //! @return Status of finding the proper transform.
  //! @retval 1: Failure.
  //! @retval 0: Success.
  friend int find_transform(const Group&, const Group&, Transform&);
  //! @brief Get the RMSD between two Groups (Conformatinos).
  //! @param Two Group's.
  //! @return RMSD.
  friend double get_rmsd(const Group&, const Group&);

  friend class Transform;
 protected:
};

Vec3d get_center_of_mass(const Group&, std::vector<double>);
double get_radius_of_gyration(const Group&);
int find_transform(const Group&, const Group&, Transform&);
double get_rmsd(const Group&, const Group&);

}

#endif

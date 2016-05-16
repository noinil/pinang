/*!
************************************************************
@file group.hpp
@brief Definition of class group.

In this file class Group is defined, which is basically same as class
Conformation, but contains more member and friend functions.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:03
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_GROUP_H
#define PINANG_GROUP_H

#include "geometry.hpp"
#include "conformation.hpp"

namespace pinang {


class Transform;

/*!
************************************************************
@brief A collection or sub-part of a conformation.

A group of coordinates, usually as a sub-group of a conformation.  Geometric or
Arithmetic claculations can be computed for different groups, such as RMSD,
superimposition, or center-of-mass distances.
************************************************************
*/
class Group : public Conformation
{
 public:
  Group();
  Group(std::vector<Vec3d>);
  virtual ~Group() {coordinates_.clear();}
  int set_conformation(std::vector<Vec3d>);

  Vec3d get_centroid() const;

  friend Vec3d get_center_of_mass(const Group&, std::vector<double>);
  friend double get_radius_of_gyration(const Group&);
  friend int find_transform(const Group&, const Group&, Transform&);
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

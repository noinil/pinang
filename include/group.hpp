// -*-c++-*-

#ifndef PINANG_GROUP_H
#define PINANG_GROUP_H

#include "geometry.hpp"
#include "conformation.hpp"

namespace pinang {

class Transform;

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
 protected:
};

Vec3d get_center_of_mass(const Group&, std::vector<double>);
double get_radius_of_gyration(const Group&);
int find_transform(const Group&, const Group&, Transform&);
double get_rmsd(const Group&, const Group&);

}

#endif

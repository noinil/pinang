// -*-c++-*-

#ifndef PINANG_GROUP_H
#define PINANG_GROUP_H

#include <vector>
#include "constants.h"
#include "conformation.h"
#include "vec3d.h"

namespace pinang {

class Group : public Conformation
{
 public:
  Group();
  Group(std::vector<Vec3d>);
  virtual ~Group() {coordinates_.clear();}

  inline int set_conformation(std::vector<Vec3d>);
  inline Vec3d get_centroid() const;
 protected:
  Vec3d centroid_;
};

Group::Group()
{
  n_atom_ = 0;
  coordinates_.clear();
  centroid_.set_coords(0.0, 0.0, 0.0);
}

Group::Group(std::vector<Vec3d> v)
{
  coordinates_ = v;
  n_atom_ = coordinates_.size();
  Vec3d com(0.0, 0.0, 0.0);
  for (int i = 0; i < n_atom_; ++i) {
    com += v[i];
  }
  com /= n_atom_;
  centroid_ = com;
}

inline Vec3d Group::get_centroid() const
{
  return centroid_;
}

int Group::set_conformation(std::vector<Vec3d> v)
{
  int m = v.size();
  if (m != n_atom_ && n_atom_ > 0)
  {
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~             PINANG :: group.h              ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cerr << " ERROR: Wrong atom number when set conformation. " << std::endl;
    return 1;
  } else {
    coordinates_ = v;
    n_atom_ = m;
    // --- compute centroid ---
    Vec3d com(0.0, 0.0, 0.0);
    for (int i = 0; i < m; ++i) {
      com += v[i];
    }
    com /= m;
    centroid_ = com;
    return 0;
  }
}

}

#endif

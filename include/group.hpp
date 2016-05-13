// -*-c++-*-

#ifndef PINANG_GROUP_H
#define PINANG_GROUP_H

#include <vector>
#include "constants.hpp"
#include "conformation.hpp"
#include "vec3d.hpp"

namespace pinang {

class Group : public Conformation
{
 public:
  Group();
  Group(std::vector<Vec3d>);
  virtual ~Group() {coordinates_.clear();}

  int set_conformation(std::vector<Vec3d>);
  Vec3d get_centroid() const;
 protected:
  Vec3d centroid_;
};

}

#endif

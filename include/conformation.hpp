#ifndef PINANG_CONFORMATION_H
#define PINANG_CONFORMATION_H

#include <vector>
#include "vec3d.hpp"

namespace pinang {

class Conformation
{
 public:
  Conformation();
  Conformation(std::vector<Vec3d>);
  virtual ~Conformation() {coordinates_.clear();}

  void reset();

  int get_size() const { return n_atom_; }

  int set_conformation(std::vector<Vec3d>);
  Vec3d& get_coor(int);

 protected:
  std::vector<Vec3d> coordinates_;
  int n_atom_;
};


}

#endif

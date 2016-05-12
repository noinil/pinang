// -*-c++-*-

#ifndef PINANG_CONFORMATION_H
#define PINANG_CONFORMATION_H

#include <vector>
#include "constants.h"
#include "vec3d.h"

namespace pinang {
class Conformation
{
 public:
  Conformation();
  Conformation(std::vector<Vec3d>);
  virtual ~Conformation() {coordinates_.clear();}

  inline void reset();

  inline int get_size() {return n_atom_;}

  int set_conformation(std::vector<Vec3d>);
  inline const Vec3d& get_coor(int) const;

 protected:
  std::vector<Vec3d> coordinates_;
  int n_atom_;
};

Conformation::Conformation()
{
  n_atom_ = 0;
  coordinates_.clear();
}

Conformation::Conformation(std::vector<Vec3d> v)
{
  coordinates_ = v;
  n_atom_ = coordinates_.size();
}

inline void Conformation::reset()
{
  n_atom_ = 0;
  coordinates_.clear();
}

int Conformation::set_conformation(std::vector<Vec3d> v)
{
  int m = v.size();
  if (m != n_atom_ && n_atom_ > 0)
  {
    std::cerr << " ERROR: Wrong atom number when set conformation. " << std::endl;
    return 1;
  } else {
    coordinates_ = v;
    n_atom_ = m;
    return 0;
  }
}

inline const Vec3d& Conformation::get_coor(int n) const
{
  if ( n >= n_atom_ || n < 0)
  {
    std::cerr << " ERROR: Atom index out of range in Conformation. " << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    return coordinates_[n];
  }
}

}

#endif

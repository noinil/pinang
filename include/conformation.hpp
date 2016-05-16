/*!
************************************************************
@file conformation.hpp
@brief Definition of class Conformation.

In this file class Conformation is defined.  Each conformation contains a list
of coordinates.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 17:57
@copyright GNU Public License V3.0
************************************************************
*/

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
  std::vector<Vec3d> coordinates_;  //!< A set of coordinate objects in a certain conformation.
  int n_atom_;                      //!< Number of coordinates in a conformation.
};


}

#endif

/*!
  @file conformation.hpp
  @brief Definition of class Conformation.

  In this file class Conformation is defined.  Each conformation contains a list
  of coordinates.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 17:57
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_CONFORMATION_H
#define PINANG_CONFORMATION_H

#include <vector>
#include "vec3d.hpp"

namespace pinang {

/*!
  @brief A certain configuration of a biomolecule.

  The class Conformation consists of a set of coordinates, ignoring the physical
  properties such as mass or charge.  Usually objects of this class is read from
  PDB or dcd files.
*/
class Conformation
{
 public:
  //! @brief Create an "empty" Conformation object.
  //! @return A Conformation object.
  Conformation();
  //! @brief Create a Conformation object based on a set of coordinates.
  //! @param A set of Vec3d coordinates.
  //! @return A Conformation object.
  Conformation(std::vector<Vec3d>);
  virtual ~Conformation() {coordinates_.clear();}

  //! @brief Reset properties of Conformation.
  void reset();

  //! @brief Get number of coordinates in Conformation.
  //! @return Number of coordinates in Conformation.
  int get_size() const { return n_atom_; }

  //! @brief Set Conformation based on a set of coordinates.
  //! @param A set of Vec3d coordinates, whose size should be same as original Conformation.
  //! @return Status of setting up the Conformation.
  //! @retval 1: Failure.
  //! @retval 0: Success.
  int set_conformation(std::vector<Vec3d>);

  //! @brief Get a coordinate from conformation.
  //! @param Serial number of the coordinate.
  //! @return Vec3d type coordinate.
  Vec3d& get_coordinate(int);

 protected:
  std::vector<Vec3d> coordinates_;  //!< A set of coordinate objects in a certain conformation.
  int n_atom_;                      //!< Number of coordinates in a conformation.
};


}

#endif

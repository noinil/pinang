/*!
  @file conformation.cpp
  @brief Define functions of class Conformation.

  Definitions of member or friend functions of class Conformation.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:38
  @copyright GNU Public License V3.0
*/

#include <fstream>
#include "conformation.hpp"

namespace pinang {

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

Conformation::Conformation(const std::string& fname)
{
  std::ifstream ifile(fname.c_str());
  if (!ifile.is_open())
  {
    std::cout << " ~       PINANG :: Conformation.cpp        ~ " << "\n";
    std::cerr << " ERROR: Cannot read file: " << fname << "\n";
    exit(EXIT_FAILURE);
  }

  n_atom_ = 0;
  Vec3d coor_tmp;
  while (ifile.good()) {
    ifile >> coor_tmp;
    if (ifile.fail())
      break;
    coordinates_.push_back(coor_tmp);
    n_atom_++;
  }
}

void Conformation::reset()
{
  n_atom_ = 0;
  coordinates_.clear();
}

int Conformation::set_conformation(std::vector<Vec3d> v)
{
  int m = v.size();
  if (m != n_atom_ && n_atom_ > 0)
  {
    std::cout << " ~             PINANG :: conformation.hpp       ~ " << "\n";
    std::cerr << " ERROR: Wrong atom number when set conformation. " << "\n";
    return 1;
  } else {
    coordinates_ = v;
    n_atom_ = m;
    return 0;
  }
}

Vec3d& Conformation::get_coordinate(int n)
{
  if (n >= n_atom_ || n < 0)
  {
    std::cout << " ~             PINANG :: conformation.hpp       ~ " << "\n";
    std::cerr << " ERROR: Atom index out of range in Conformation. " << "\n";
    exit(EXIT_SUCCESS);
  } else {
    return coordinates_[n];
  }
}

}  // pinang

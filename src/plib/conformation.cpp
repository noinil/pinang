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
    std::cout << " ~             PINANG :: conformation.hpp       ~ " << std::endl;
    std::cerr << " ERROR: Wrong atom number when set conformation. " << std::endl;
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
    std::cout << " ~             PINANG :: conformation.hpp       ~ " << std::endl;
    std::cerr << " ERROR: Atom index out of range in Conformation. " << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    return coordinates_[n];
  }
}

}  // pinang

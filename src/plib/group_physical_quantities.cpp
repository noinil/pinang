#include "group.hpp"

namespace pinang {

Vec3d get_center_of_mass(const Group& grp, std::vector<double> masses)
{
  int m1 = grp.n_atom_;
  int m2 = masses.size();
  if (m1 != m2) {
    std::cout << " ERROR: inconsistent number of atoms when calculating Rg ";
    exit(EXIT_SUCCESS);
  }

  return Vec3d(0, 0, 0);
}


}  // pinang

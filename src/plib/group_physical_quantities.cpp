#include <cmath>
#include "group.hpp"

namespace pinang {

Vec3d get_center_of_mass(const Group& grp, std::vector<double> masses)
{
  Vec3d com;
  double M = 0;  // sum of masses;

  int m1 = grp.n_atom_;
  int m2 = masses.size();
  if (m1 != m2) {
    std::cout << " ERROR: inconsistent number of atoms when calculating Rg! (group.hpp)\n";
    exit(EXIT_SUCCESS);
  }
  if (m2 == 0) {
    std::cout << " ERROR: sum of masses == 0!  (group.hpp)\n";
    exit(EXIT_SUCCESS);
  }

  for (int i = 0; i < m1; ++i) {
    com += masses[i] * grp.coordinates_[i];
    M += masses[i];
  }
  com /= M;

  return com;
}

double get_radius_of_gyration(const Group& grp)
{
  double rg = 0;
  Vec3d ctr = grp.centroid_;
  int m = grp.n_atom_;

  Vec3d vtmp;
  double dtmp;
  double dtmp_sum = 0;

  for (int i = 0; i < m; ++i) {
    vtmp = grp.coordinates_[i] - ctr;
    dtmp = vtmp.squared_norm();
    dtmp_sum += dtmp;
  }
  rg = sqrt(dtmp_sum / m);

  return rg;
}


}  // pinang

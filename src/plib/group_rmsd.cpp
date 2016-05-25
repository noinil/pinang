/*!
  @file group_rmsd.cpp
  @brief Define function of RMSD calculation.

  Definition of function of superimposition transform and RMSD calculation.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:43
  @copyright GNU Public License V3.0
*/


// Based on https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp;

#include <cmath>
#include <Eigen/Geometry>
#include "group.hpp"

namespace pinang {

int find_transform(const Group& grp1, const Group& grp2, Transform& t)
{
  // Simple check...
  int m1 = grp1.n_atom_;
  int m2 = grp2.n_atom_;
  if (m1 != m2) {
    std::cout << " ERROR: inconsistent number of atoms when calculating Rg! (group_rmsd.cpp)\n";
    exit(EXIT_SUCCESS);
  }

  // Make copies of input groups...
  Group g1 = grp1;
  Group g2 = grp2;

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  double dist_in = 0, dist_out = 0;
  Vec3d vtmp;
  for (int i = 0; i < m1-1; ++i) {
    vtmp = g1.coordinates_[i + 1] - g1.coordinates_[i];
    dist_in += vtmp.norm();
    vtmp = g2.coordinates_[i + 1] - g2.coordinates_[i];
    dist_out += vtmp.norm();
  }
  if (dist_in <= 0 || dist_out <= 0)
    return 1;  // 1 means failure of superimposition;
  double scale = dist_out/dist_in;
  for (int i = 0; i < m1; ++i) {
    g2.coordinates_[i] /= scale;
  }
  Vec3d pctr1 = g1.get_centroid();
  Vec3d pctr2 = g2.get_centroid();
  Eigen::Vector3d in_ctr(pctr1[0], pctr1[1], pctr1[2]);
  Eigen::Vector3d out_ctr(pctr2[0], pctr2[1], pctr2[2]);

  // Find the centroids then shift to the origin
  for (int i = 0; i < m1; ++i) {
    g1.coordinates_[i] -= pctr1;
    g2.coordinates_[i] -= pctr2;
  }

  // pinang::Group ---> Eigen::Matrix
  int col = 0;
  int row = 0;
  Eigen::Matrix3Xd in(3, m1);
  Eigen::Matrix3Xd out(3, m1);
  for (col = 0; col < m1; ++col) {
    in(0, col) = g1.coordinates_[col][0];
    in(1, col) = g1.coordinates_[col][1];
    in(2, col) = g1.coordinates_[col][2];
    out(0, col) = g2.coordinates_[col][0];
    out(1, col) = g2.coordinates_[col][1];
    out(2, col) = g2.coordinates_[col][2];
  }

  // ==================== KEY!!! ====================
  // SVD
  Eigen::MatrixXd Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  I(2, 2) = d;
  Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

  // prepare for the transformation vector
  Eigen::Quaterniond eigen_q(R);  // rotation part!
  Eigen::Vector3d eigen_t;
  eigen_t = out_ctr * scale - R * in_ctr;  // translation part!

  // ==================== SAVE to TRANSFORM VECTOR ====================
  Quaternion transform_linear_q(eigen_q.w(), eigen_q.x(), eigen_q.y(), eigen_q.z());
  Vec3d transform_translation(eigen_t(0), eigen_t(1), eigen_t(2));
  t.set_rotation(transform_linear_q);
  t.set_translation(transform_translation);

  return 1;
}

double get_rmsd(const Group& grp1, const Group& grp2)
{
  double rmsd = 0;
  double d = 0;
  Vec3d vtmp;

  // Simple check...
  int m1 = grp1.n_atom_;
  int m2 = grp2.n_atom_;
  if (m1 != m2) {
    std::cout << " ERROR: inconsistent number of atoms when calculating Rg! (group_rmsd.cpp)\n";
    exit(EXIT_SUCCESS);
  }

  // rmsd calculation loop
  for (int i = 0; i < m1; ++i) {
    vtmp = grp1.coordinates_[i] - grp2.coordinates_[i];
    d += vtmp.squared_norm();
  }

  rmsd = sqrt(d / m1);
  return rmsd;
}

}  // pinang

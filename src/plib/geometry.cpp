#include <cmath>
#include "geometry.hpp"

namespace pinang {

int Quaternion::normalize()
{
  double d = w_ * w_ + x_ * x_ + y_ * y_ + z_ * z_;
  d = sqrt(d);
  if (d <= 0) {
    std::cout << " ~             PINANG :: geometry.cpp       ~ " << std::endl;
    std::cerr << " WARNING: wrong quaternion (equals to 0, 0, 0, 0). " << std::endl;
    return 1;
  }
  w_ /= d;
  x_ /= d;
  y_ /= d;
  z_ /= d;
  return 0;
}

Transform::Transform()
{
  Quaternion q;
  Vec3d v(0.0, 0.0, 0.0);
  rotation_ = q;
  translation_ = v;
  rotv1_ = Vec3d(1.0, 0.0, 0.0);
  rotv2_ = Vec3d(0.0, 1.0, 0.0);
  rotv3_ = Vec3d(0.0, 0.0, 1.0);
}

Transform::Transform(const Quaternion& q, const Vec3d& v)
{
  rotation_ = q;
  translation_ = v;
  rotation_.normalize();
  double w = rotation_.w();
  double x = rotation_.x();
  double y = rotation_.y();
  double z = rotation_.z();
  double x2 = 2 * x * x;
  double y2 = 2 * y * y;
  double z2 = 2 * z * z;
  double xw = 2 * x * w;
  double yw = 2 * y * w;
  double zw = 2 * z * w;
  double xy = 2 * x * y;
  double yz = 2 * z * y;
  double zx = 2 * x * z;
  rotv1_ = Vec3d(1 - y2 - z2, xy - zw, zx + yw);
  rotv2_ = Vec3d(xy + zw, 1 - x2 - z2, yz - xw);
  rotv3_ = Vec3d(zx - yw, yz + xw, 1 - x2 - y2);
}

void Transform::set_rotation(const Quaternion& q)
{
  rotation_ = q;
  rotation_.normalize();
  double w = rotation_.w();
  double x = rotation_.x();
  double y = rotation_.y();
  double z = rotation_.z();
  double x2 = 2 * x * x;
  double y2 = 2 * y * y;
  double z2 = 2 * z * z;
  double xw = 2 * x * w;
  double yw = 2 * y * w;
  double zw = 2 * z * w;
  double xy = 2 * x * y;
  double yz = 2 * z * y;
  double zx = 2 * x * z;
  rotv1_ = Vec3d(1 - y2 - z2, xy - zw, zx + yw);
  rotv2_ = Vec3d(xy + zw, 1 - x2 - z2, yz - xw);
  rotv3_ = Vec3d(zx - yw, yz + xw, 1 - x2 - y2);
}

Vec3d Transform::apply(const Vec3d& v)
{
  Vec3d v_out;
  double x1 = rotv1_ * v;
  double y1 = rotv2_ * v;
  double z1 = rotv3_ * v;
  v_out.set_coords(x1, y1, z1);
  v_out += translation_;

  return v_out;
}

Group Transform::apply(const Group& g)
{
  std::vector<Vec3d> vv;
  Vec3d vtmp;
  double x1, y1, z1;
  int s = g.n_atom_;
  for (int i = 0; i < s; ++i) {
    vtmp = g.coordinates_[i];
    x1 = rotv1_ * vtmp;
    y1 = rotv2_ * vtmp;
    z1 = rotv3_ * vtmp;
    vtmp.set_coords(x1, y1, z1);
    vtmp += translation_;
    vv.push_back(vtmp);
  }
  Group gtmp(vv);

  return gtmp;
}


}  // pinang

#include <sstream>
#include "particle.hpp"

namespace pinang {

void Particle::set_atom_name(const std::string& s)
{
  size_t sz = 4;
  atom_name_ = s;

  if (atom_name_.size() < sz)
  {
    atom_name_.resize(sz, ' ');
  }
}

void Particle::set_resid_name(const std::string& s)
{
  size_t sz = 3;
  resid_name_ = s;

  if (resid_name_.size() < sz)
  {
    resid_name_.resize(sz, ' ');
  }
}

// Particle::Particle ==============================================================
Particle::Particle()
{
  atom_name_ = "";
  resid_name_ = "";
  resid_index_ = 0;
  charge_ = 0;
  mass_ = 0;
}

void Particle::reset()
{
  atom_name_ = "";
  resid_name_ = "";
  resid_index_ = 0;
  charge_ = 0;
  mass_ = 0;
}


std::istream& operator>>(std::istream& i, Particle& p)
{
  std::istringstream tmp_sstr;
  std::string top_line;

  std::string tmp_s;
  int tmp_i;
  double tmp_d;

  std::getline(i, top_line);
  tmp_sstr.str(top_line);

  tmp_sstr >> tmp_i;

  tmp_sstr >> tmp_i;
  p.resid_index_ = tmp_i;

  tmp_sstr >> tmp_s;
  p.set_resid_name(tmp_s);

  tmp_sstr >> tmp_s;
  p.set_atom_name(tmp_s);

  tmp_sstr >> tmp_d;
  p.mass_ = tmp_d;

  tmp_sstr >> tmp_d;
  p.charge_ = tmp_d;

  if (!i) return i;
  return i;
}

}  // pinang

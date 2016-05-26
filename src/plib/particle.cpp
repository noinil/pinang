/*!
  @file particle.cpp
  @brief Define functions of class Particle.

  Definitions of member or friend functions of class Particle.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:44
  @copyright GNU Public License V3.0
*/

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

void Particle::set_residue_name(const std::string& s)
{
  size_t sz = 3;
  residue_name_ = s;

  if (residue_name_.size() < sz)
  {
    residue_name_.resize(sz, ' ');
  }
}

// Particle::Particle ==============================================================
Particle::Particle()
{
  atom_name_ = "";
  residue_name_ = "";
  residue_serial_ = 0;
  charge_ = 0;
  mass_ = 0;
}

void Particle::reset()
{
  atom_name_ = "";
  residue_name_ = "";
  residue_serial_ = 0;
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
  char tmp_c;

  std::getline(i, top_line);
  tmp_sstr.str(top_line);

  tmp_sstr >> tmp_i;

  tmp_sstr >> tmp_c;
  p.chain_ID_ = tmp_c;

  tmp_sstr >> tmp_i;
  p.residue_serial_ = tmp_i;

  tmp_sstr >> tmp_s;
  p.set_residue_name(tmp_s);

  tmp_sstr >> tmp_s;
  p.set_atom_name(tmp_s);

  tmp_sstr >> tmp_s;

  tmp_sstr >> tmp_d;
  p.charge_ = tmp_d;

  tmp_sstr >> tmp_d;
  p.mass_ = tmp_d;

  if (!i) return i;
  return i;
}

}  // pinang

/*!
  @file residue.cpp
  @brief Define functions of class Residue.

  Definitions of member or friend functions of class Residue.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:45
  @copyright GNU Public License V3.0
*/


#include "residue.hpp"

namespace pinang {

void Residue::set_residue_name(const std::string& s)
{
  size_t sz = 3;
  residue_name_ = s;
  if (residue_name_.size() < sz)
  {
    residue_name_.resize(sz, ' ');
  }
}
void Residue::set_residue_by_name(const std::string& s)
{
  PhysicalProperty p;
  size_t sz = 3;
  residue_name_ = s;
  if (residue_name_.size() < sz)
  {
    residue_name_.resize(sz, ' ');
  }
  short_name_ = p.get_short_name(s);
  chain_type_ = p.get_chain_type(s);
  mass_ = p.get_mass(s);
  charge_ = p.get_charge(s);
}

Atom& Residue::get_atom(int n)
{
  if (v_atoms_.empty())
  {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: No Atoms in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  if (n >= int(v_atoms_.size()))
  {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: Atom index out of range in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return v_atoms_[n];
}


int Residue::add_atom(const Atom& a)
{
  if (a.get_residue_serial() != residue_serial_)
  {
    return 1;
  }
  for (const Atom& b : v_atoms_) {
    if (a.get_atom_name() == b.get_atom_name())
      return 0;
  }  // in case of NMR uncertain multi atoms
  v_atoms_.push_back(a);

  std::string an = a.get_atom_name();
  if (an == "CA  ")
  {
    cg_C_alpha_ = a;
  }
  if (an == "CB  ")
  {
    cg_C_beta_ = a;
  }
  if (an == "C3' " || an == "S   " || an == "DS  ")
  {
    cg_S_ = a;
  }
  if (an == "P   " || an == "DP  ")
  {
    cg_P_ = a;
  }
  if (an == "N1  " || an == "B   " || an == "DB  ")
  {
    cg_B_ = a;
  }
  if (a.get_record_name() == "HETATM" && a.get_element() != "H")
  {
    cg_C_alpha_ = a;
  }
  ++n_atom_;
  return 0;
}

int Residue::delete_atom(const int i)
{
  if (i >= n_atom_){
    return 1;
  }
  v_atoms_.erase(v_atoms_.begin() + i);
  n_atom_--;
  return 0;
}



Atom& Residue::get_cg_C_alpha()
{
  if (cg_C_alpha_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: C_alpha not set in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return cg_C_alpha_;
}

Atom& Residue::get_cg_C_beta()
{
  if (cg_C_beta_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: C_beta not set in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return cg_C_beta_;
}

Atom& Residue::get_cg_P()
{
  if (cg_P_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: CG Phosphate not set in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return cg_P_;
}

Atom& Residue::get_cg_S()
{
  if (cg_S_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: CG Sugar not set in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return cg_S_;
}

Atom& Residue::get_cg_B()
{
  if (cg_B_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
    std::cerr << "ERROR: CG Base not set in Residue: "
              << residue_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return cg_B_;
}

void Residue::set_cg_C_alpha()
{
  for (const Atom& b : v_atoms_) {
    if (b.get_atom_name() == "CA  ")
      cg_C_alpha_ = b;
    break;
  }
}

void Residue::set_cg_C_beta()
{
  for (const Atom& b : v_atoms_) {
    if (b.get_atom_name() == "CB  ")
      cg_C_beta_ = b;
    break;
  }
}

void Residue::self_check() const
{
  for (const Atom& a : v_atoms_) {
    if (a.get_chain_ID() != chain_ID_ || a.get_residue_serial() != residue_serial_
        || a.get_residue_name() != residue_name_)
    {
      std::cout << " ~               PINANG :: residues.hpp         ~ " << "\n";
      std::cerr << "ERROR: Inconsistent chain ID or residue index or residue type in Residue "
                << residue_serial_ << "\n";
      exit(EXIT_SUCCESS);
    }
  }
}


Residue::Residue()
{
  residue_name_ = "";
  short_name_ = "0";
  chain_ID_ = -1;
  residue_serial_ = -1;
  v_atoms_.clear();
  n_atom_ = 0;
  charge_ = 0.0;
  mass_ = 100.0;
  terminus_flag_ = 0;

  cg_C_alpha_.reset();
  cg_C_beta_.reset();
  cg_P_.reset();
  cg_S_.reset();
  cg_B_.reset();
  chain_type_ = none;
}

void Residue::reset()
{
  residue_name_ = "";
  short_name_ = "0";
  chain_ID_ = -1;
  residue_serial_ = -1;
  v_atoms_.clear();
  n_atom_ = 0;
  charge_ = 0.0;
  mass_ = 100.0;
  terminus_flag_ = 0;

  chain_type_ = none;
  cg_P_.reset();
  cg_S_.reset();
  cg_B_.reset();
  cg_C_alpha_.reset();
  cg_C_beta_.reset();
}


std::ostream& operator<<(std::ostream& o, Residue& r)
{
  int i = 0;
  int s = r.n_atom_;
  for (i = 0; i < s; ++i) {
    o << r.v_atoms_[i];
  }
  return o;
}

double residue_min_distance(const Residue& r1, const Residue& r2)
{
  int i, j;
  double d = atom_distance(r1.v_atoms_[0], r2.v_atoms_[0]);  // min_distance;
  double f = 0.0;           // tmp distance;
  int s1 = r1.n_atom_;
  int s2 = r2.n_atom_;
  for (i = 0; i < s1; ++i) {
    if (r1.v_atoms_[i].get_element() == "H")
      continue;
    for (j = 0; j < s2; ++j) {
      if (r2.v_atoms_[j].get_element() == "H")
        continue;
      f = atom_distance(r1.v_atoms_[i], r2.v_atoms_[j]);
      if (d > f)
        d = f;
    }
  }
  return d;
}

double residue_ca_distance(const Residue& r1, const Residue& r2)
{
  double d = -1;           // distance;
  d = atom_distance(r1.cg_C_alpha_, r2.cg_C_alpha_);
  return d;
}

}  // pinang

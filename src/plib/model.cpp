/*!
  @file model.cpp
  @brief Define functions of class Model.

  Definitions of member or friend functions of class Model.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:43
  @copyright GNU Public License V3.0
*/

#include <iomanip>
#include "model.hpp"

namespace pinang {

Chain& Model::get_chain(unsigned int n)
{
  if (v_chains_.empty())
  {
    std::cout << " ~             PINANG :: model.hpp              ~ " << "\n";
    std::cerr << "ERROR: No Chains found in Model: "
              << model_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  if (n >= v_chains_.size())
  {
    std::cout << " ~             PINANG :: model.hpp              ~ " << "\n";
    std::cerr << "ERROR: Chain number out of range in Model: "
              << model_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return v_chains_[n];
}

void Model::add_chain(Chain& c)
{
  c.self_check();
  v_chains_.push_back(c);
  ++n_chain_;
}


void Model::output_sequence(int n) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    v_chains_[i].output_sequence(n);
}

void Model::output_sequence_fasta(std::ostream & f_fasta, std::string s) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    v_chains_[i].output_sequence_fasta(f_fasta, s);
}

// Model ===================================================================
Model::Model()
{
  model_serial_ = 0;
  v_chains_.clear();
  n_chain_ = 0;
}

void Model::reset()
{
  model_serial_ = 0;
  v_chains_.clear();
  n_chain_ = 0;
}

void Model::output_cg_crd(std::ostream& o)
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_cg_crd(o);
  }
}

void Model::output_cg_pdb(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_cg_pdb(o, n, m);
  }
  o << "ENDMDL" << std::endl;
}

void Model::output_top_mass(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 3 - 1;
    else
      n += v_chains_[i].get_size();
  }

  o << std::setw(8) << n << " !NATOM \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none) {
      continue;
    }
    v_chains_[i].output_top_mass(o, n, m);
  }
  o << std::endl;
}

void Model::output_top_bond(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none || ct == ion)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 3 - 2;
    else
      n += v_chains_[i].get_size() - 1;
  }

  o << std::setw(8) << n << " !NBOND: bonds \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_top_bond(o, n, m);
  }
  o << " \n" << std::endl;
}

void Model::output_top_angle(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none || ct == ion)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 4 - 5;
    else
      n += v_chains_[i].get_size() - 2;
  }

  o << std::setw(8) << n << " !NTHETA: angles \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_top_angle(o, n, m);
  }
  o << " \n" << std::endl;
}

void Model::output_top_dihedral(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none || ct == ion)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 2 - 4;
    else
      n += v_chains_[i].get_size()-3;
  }

  o << std::setw(8) << n << " !NPHI: dihedrals \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_top_dihedral(o, n, m);
  }
  o << "\n" << std::endl;
}

void Model::output_ffparm_bond(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 3 - 2;
    else
      n += v_chains_[i].get_size() - 1;
  }

  o << "[ bonds ]" << std::setw(8) << n << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj"
    << std::setw(17) << "r0" << std::setw(9) << "K_b" << "\n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_ffparm_bond(o, n);
  }
  o << " \n" << std::endl;
}

void Model::output_ffparm_angle(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 4 - 5;
    else
      n += v_chains_[i].get_size()-2;
  }

  o << "[ angles ]" << std::setw(8) << n << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj" << std::setw(9) << "pk"
    << std::setw(13) << "theta_0" << std::setw(9) << "K_a" << "\n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_ffparm_angle(o, n);
  }
  o << "\n" << std::endl;
}

void Model::output_ffparm_dihedral(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 2 - 4;
    else
      n += v_chains_[i].get_size()-3;
  }

  o << "[ dihedrals ]" << std::setw(8) << n << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj"
    << std::setw(9) << "pk" << std::setw(9) << "pl"
    << std::setw(13) << "phi_0" << std::setw(9) << "K_d_1"
    << std::setw(9) << "K_d_3" << "\n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_ffparm_dihedral(o, n);
  }
  o << "\n" << std::endl;
}

std::ostream& operator<<(std::ostream& o, Model& m)
{
  o << "MODEL "
    << std::setw(8) << m.model_serial_
    << "\n";
  int i = 0;
  int s = m.n_chain_;
  for (i = 0; i < s; ++i) {
    o << m.v_chains_[i];
  }
  o << "ENDMDL" << "\n";
  return o;
}

}  // pinang

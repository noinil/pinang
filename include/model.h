// -*-c++-*-

#ifndef PINANG_MODEL_H_
#define PINANG_MODEL_H_

#include <iostream>
#include "chain.h"

namespace pinang {

class Model
{
 public:
  Model();
  virtual ~Model() {chains_.clear();}

  inline void reset();

  inline int get_model_ID() const;
  inline void set_model_ID(int n);

  inline Chain& get_chain(unsigned int n);
  inline void add_chain(Chain& c);

  inline void print_sequence(int n) const;
  inline void output_fasta(std::ostream & f_fasta, std::string s) const;

  inline int get_model_size() const;

  void output_cg_pos(std::ostream& o);

  void output_top_mass(std::ostream& o);
  void output_top_bond(std::ostream& o);
  void output_top_angle(std::ostream& o);
  void output_top_dihedral(std::ostream& o);

  void output_top_nonbonded(std::ostream& o);

 protected:
  int model_ID_;
  std::vector<Chain> chains_;
  int n_chain_;
};

inline int Model::get_model_ID() const
{
  return model_ID_;
}
inline void Model::set_model_ID(int n)
{
  model_ID_ = n;
}


inline Chain& Model::get_chain(unsigned int n)
{
  if (chains_.empty())
  {
    std::cerr << "ERROR: No Chains found in Model: "
              << model_ID_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= chains_.size())
  {
    std::cerr << "ERROR: Chain number out of range in Model: "
              << model_ID_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return chains_[n];
}

inline void Model::add_chain(Chain& c)
{
  c.self_check();
  chains_.push_back(c);
  n_chain_++;
}


inline void Model::print_sequence(int n) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    chains_[i].pr_seq(n);
}

inline void Model::output_fasta(std::ostream & f_fasta, std::string s) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    chains_[i].output_fasta(f_fasta, s);
}


inline int Model::get_model_size() const
{
  return n_chain_;
}

// Model ===================================================================
inline Model::Model()
{
  model_ID_ = 0;
  chains_.clear();
  n_chain_ = 0;
}

inline void Model::reset()
{
  model_ID_ = 0;
  chains_.clear();
  n_chain_ = 0;
}


void Model::output_cg_pos(std::ostream& o)
{
  int i = 0;
  int n = 0;
  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    chains_[i].output_cg_pos(o, n);
  }
}

void Model::output_top_mass(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += chains_[i].get_chain_length() * 3 - 1;
    else
      n += chains_[i].get_chain_length();
  }

  o << "[ particles ]"
    << std::setw(6) << n
    << std::endl;
  o << "# "
    << std::setw(9) << "index"
    << std::setw(8) << "resid"
    << std::setw(8) << "resname"
    << std::setw(8) << "atom"
    << std::setw(10) << "mass"
    << std::setw(8) << "charge"
    << std::endl;

  n = 0;
  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    chains_[i].output_top_mass(o, n);
  }
  o << std::endl;
}

void Model::output_top_bond(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += chains_[i].get_chain_length() * 3 - 2;
    else
      n += chains_[i].get_chain_length() - 1;
  }

  o << "[ bonds ]"
    << std::setw(6) << n
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(6) << "pj"
    << std::setw(10) << "r0"
    << std::setw(8) << "K_b"
    << std::endl;

  n = 0;
  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    chains_[i].output_top_bond(o, n);
  }
  o << std::endl;
}

void Model::output_top_angle(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += chains_[i].get_chain_length() * 4 - 5;
    else
      n += chains_[i].get_chain_length()-2;
  }

  o << "[ angles ]"
    << std::setw(6) << n
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(6) << "pj"
    << std::setw(6) << "pk"
    << std::setw(12) << "theta_0"
    << std::setw(8) << "K_a"
    << std::endl;

  n = 0;
  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    chains_[i].output_top_angle(o, n);
  }
  o << std::endl;
}

void Model::output_top_dihedral(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += chains_[i].get_chain_length() * 2 - 4;
    else
      n += chains_[i].get_chain_length()-3;
  }

  o << "[ dihedrals ]"
    << std::setw(6) << n
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(6) << "pj"
    << std::setw(6) << "pk"
    << std::setw(6) << "pl"
    << std::setw(12) << "phi_0"
    << std::setw(8) << "K_d_1"
    << std::setw(8) << "K_d_3"
    << std::endl;

  n = 0;
  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    chains_[i].output_top_dihedral(o, n);
  }
  o << std::endl;
}

void Model::output_top_nonbonded(std::ostream& o)
{
  int i = 0;
  Chain c0;
  Chain c_tmp;
  Residue r_tmp;

  for (i = 0; i < n_chain_; i++) {
    ChainType ct = chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    if (ct == DNA || ct == RNA || ct == na)
    {
      c_tmp.reset();

      Atom S = chains_[i].get_residue(0).get_S();
      Atom B = chains_[i].get_residue(0).get_B();
      r_tmp.reset();
      r_tmp.set_residue_by_name(S.get_resid_name());
      r_tmp.set_chain_ID(S.get_chain_ID());
      r_tmp.set_resid_index(S.get_resid_index());
      r_tmp.add_atom(S);
      c_tmp.add_residue(r_tmp);

      r_tmp.reset();
      r_tmp.set_residue_by_name(B.get_resid_name());
      r_tmp.set_chain_ID(B.get_chain_ID());
      r_tmp.set_resid_index(B.get_resid_index());
      r_tmp.add_atom(B);
      c_tmp.add_residue(r_tmp);


      for (int j = 1; j < chains_[i].get_chain_length(); j++) {
        Atom P = chains_[i].get_residue(j).get_P();
        Atom S = chains_[i].get_residue(j).get_S();
        Atom B = chains_[i].get_residue(j).get_B();
        r_tmp.reset();
        r_tmp.set_residue_by_name(P.get_resid_name());
        r_tmp.set_chain_ID(P.get_chain_ID());
        r_tmp.set_resid_index(P.get_resid_index());
        r_tmp.add_atom(P);
        c_tmp.add_residue(r_tmp);

        r_tmp.reset();
        r_tmp.set_residue_by_name(S.get_resid_name());
        r_tmp.set_chain_ID(S.get_chain_ID());
        r_tmp.set_resid_index(S.get_resid_index());
        r_tmp.add_atom(S);
        c_tmp.add_residue(r_tmp);

        r_tmp.reset();
        r_tmp.set_residue_by_name(B.get_resid_name());
        r_tmp.set_chain_ID(B.get_chain_ID());
        r_tmp.set_resid_index(B.get_resid_index());
        r_tmp.add_atom(B);
        c_tmp.add_residue(r_tmp);
      }
      c_tmp.set_chain_type(ct);
      c0 = c0 + c_tmp;
      continue;
    }
    c0 = c0 + chains_[i];
  }

  o << "[ native ]"
    << std::setw(6) << c0.get_native_contact_number()
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(6) << "pj"
    << std::setw(8) << "eps"
    << std::setw(10) << "sigma"
    << std::endl;

  c0.output_top_native(o);
  o << std::endl;
}


inline std::ostream& operator<<(std::ostream& o, Model& m)
{
  o << "MODEL "
    << std::setw(8) << m.get_model_ID()
    << std::endl;
  int i = 0;
  for (i = 0; i < m.get_model_size(); i++) {
    o << m.get_chain(i) ;
  }
  o << "ENDMDL" << std::endl;
  return o;
}

}

#endif

#include <iomanip>
#include "model.hpp"

namespace pinang {

Chain& Model::get_chain(unsigned int n)
{
  if (chains_.empty())
  {
    std::cout << " ~             PINANG :: model.hpp              ~ " << std::endl;
    std::cerr << "ERROR: No Chains found in Model: "
              << model_ID_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= chains_.size())
  {
    std::cout << " ~             PINANG :: model.hpp              ~ " << std::endl;
    std::cerr << "ERROR: Chain number out of range in Model: "
              << model_ID_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return chains_[n];
}

void Model::add_chain(Chain& c)
{
  c.self_check();
  chains_.push_back(c);
  n_chain_++;
}


void Model::print_sequence(int n) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    chains_[i].pr_seq(n);
}

void Model::output_fasta(std::ostream & f_fasta, std::string s) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    chains_[i].output_fasta(f_fasta, s);
}

// Model ===================================================================
Model::Model()
{
  model_ID_ = 0;
  chains_.clear();
  n_chain_ = 0;
}

void Model::reset()
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
    << std::setw(8) << n
    << std::endl;
  o << "# "
    << std::setw(9) << "index"
    << std::setw(9) << "resid"
    << std::setw(10) << "resname"
    << std::setw(10) << "atom"
    << std::setw(17) << "mass"
    << std::setw(13) << "charge"
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
    << std::setw(8) << n
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(9) << "pj"
    << std::setw(17) << "r0"
    << std::setw(9) << "K_b"
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
    << std::setw(8) << n
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(9) << "pj"
    << std::setw(9) << "pk"
    << std::setw(13) << "theta_0"
    << std::setw(9) << "K_a"
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
    << std::setw(8) << n
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(9) << "pj"
    << std::setw(9) << "pk"
    << std::setw(9) << "pl"
    << std::setw(13) << "phi_0"
    << std::setw(9) << "K_d_1"
    << std::setw(9) << "K_d_3"
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
    << std::setw(8) << c0.get_native_contact_number()
    << std::endl;
  o << "# "
    << std::setw(6) << "pi"
    << std::setw(9) << "pj"
    << std::setw(17) << "sigma"
    << std::setw(13) << "eps"
    << std::endl;

  c0.output_top_native(o);
  o << std::endl;
}


std::ostream& operator<<(std::ostream& o, Model& m)
{
  o << "MODEL "
    << std::setw(8) << m.model_ID_
    << std::endl;
  int i = 0;
  int s = m.n_chain_;
  for (i = 0; i < s; i++) {
    o << m.chains_[i];
  }
  o << "ENDMDL" << std::endl;
  return o;
}

}  // pinang

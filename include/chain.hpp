#ifndef PINANG_CHAIN_H_
#define PINANG_CHAIN_H_

#include <iostream>
#include "residue.hpp"

namespace pinang {

class Chain
{
 public:
  Chain();
  virtual ~Chain() {residues_.clear();}

  void reset();

  char get_chain_ID() const { return chain_ID_; }
  void set_chain_ID(char a) { chain_ID_ = a; }

  ChainType get_chain_type() const { return chain_type_; }
  void set_chain_type(ChainType a) { chain_type_ = a; }

  Residue& get_residue(int);
  int add_residue(const Residue&);

  int get_chain_length() const { return n_residue_; }

  void pr_seq(int) const;
  void output_fasta(std::ostream&, std::string) const;
  void self_check();

  void output_cg_pos(std::ostream&, int&);
  void output_top_mass(std::ostream&, int&);
  void output_top_bond(std::ostream&, int&);
  void output_top_angle(std::ostream&, int&);
  void output_top_dihedral(std::ostream&, int&);
  void output_top_native(std::ostream&);
  int get_native_contact_number();

  friend Chain operator+(const Chain&, const Chain&);

  friend std::ostream& operator<<(std::ostream&, Chain&);

 protected:
  char chain_ID_;
  ChainType chain_type_;
  int n_residue_;
  std::vector<Residue> residues_;
};

}

#endif

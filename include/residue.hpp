#ifndef PINANG_RESIDUE_H_
#define PINANG_RESIDUE_H_

#include <iostream>
#include <vector>

#include "atom.hpp"
#include "constants.hpp"

namespace pinang {

class Residue
{
 public:
  Residue();
  virtual ~Residue() {atoms_.clear();}

  void reset();

  std::string get_resid_name() const { return resid_name_; }
  std::string get_short_name() const { return short_name_; }
  void set_resid_name(const std::string&);
  void set_residue_by_name(const std::string&);

  char get_chain_ID() const { return chain_ID_; }
  void set_chain_ID(char a) { chain_ID_ = a; }

  ChainType get_chain_type() const { return chain_type_; }
  void set_chain_type(ChainType a) { chain_type_ = a; }

  int get_resid_index() const { return resid_index_; }
  void set_resid_index(int i) { resid_index_ = i; }

  int get_term_flag() const { return term_flag_; }
  // 5: 5';   3: 3';   0: not terminus;
  // -1: N;   1: C;
  void set_term_flag(int i) { term_flag_ = i; }

  double get_resid_charge() const { return charge_; }
  void set_resid_charge(double c) { charge_ = c; }

  double get_resid_mass() const { return mass_; }
  void set_resid_mass(double m) { mass_ = m; }

  void self_check() const;

  Atom& get_atom(int);
  int add_atom(const Atom&);
  int delete_atom(const int);

  int get_residue_size() const { return n_atom_; }

  Atom& get_C_alpha();
  Atom& get_C_beta();
  Atom& get_P();
  Atom& get_S();
  Atom& get_B();
  void set_C_alpha();
  void set_C_beta();
  void set_P(const Atom& a) { P_ = a; }
  void set_S(const Atom& a) { S_ = a; }
  void set_B(const Atom& a) { B_ = a; }

  friend std::ostream& operator<<(std::ostream&, Residue&);
  friend double resid_min_distance(const Residue&, const Residue&);
  friend double resid_ca_distance(const Residue&, const Residue&);
 protected:
  std::string resid_name_;
  std::string short_name_;
  char chain_ID_;
  int resid_index_;
  std::vector<Atom> atoms_;
  int n_atom_;
  double charge_;
  double mass_;

  Atom C_alpha_;
  Atom C_beta_;
  Atom P_;
  Atom S_;
  Atom B_;
  ChainType chain_type_;

  int term_flag_;  // 5: 5'; 3: 3'; -1: N;  1: C;  0: not terminus;
};

double resid_min_distance(const Residue&, const Residue&);
double resid_ca_distance(const Residue&, const Residue&);
}

#endif

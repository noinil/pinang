// -*-c++-*-

#ifndef PINANG_RESIDUE_H_
#define PINANG_RESIDUE_H_

#include <iostream>
#include <vector>
#include <cstdlib>

#include "vec3d.h"
#include "atom.h"

namespace pinang {

class Residue
{
 public:
  Residue();
  virtual ~Residue() {atoms_.clear();}

  inline void reset();

  inline std::string get_resid_name() const;
  inline std::string get_short_name() const;
  void set_resid_name(const std::string& s);
  void set_residue_by_name(const std::string& s);

  inline char get_chain_ID() const;
  inline void set_chain_ID(char a);

  inline ChainType get_chain_type() const;
  inline void set_chain_type(ChainType ct);

  inline int get_resid_index() const;
  inline void set_resid_index(int i);

  inline int get_term_flag() const;
  inline void set_term_flag(int i);

  inline double get_resid_charge() const;
  inline void set_resid_charge(double c);

  inline double get_resid_mass() const;
  inline void set_resid_mass(double c);

  inline void self_check() const;

  inline Atom& get_atom(int n);
  inline int add_atom(const Atom& a);
  inline int delete_atom(const int i);

  inline int get_residue_size() const;

  inline Atom& get_C_alpha();
  inline Atom& get_C_beta();
  inline Atom& get_P();
  inline Atom& get_S();
  inline Atom& get_B();
  inline void set_C_alpha();
  inline void set_C_beta();
  inline void set_P(const Atom& a);
  inline void set_S(const Atom& a);
  inline void set_B(const Atom& a);

  // void set_cg_na();

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


inline std::string Residue::get_resid_name() const
{
  return resid_name_;
}
inline std::string Residue::get_short_name() const
{
  return short_name_;
}


void Residue::set_resid_name(const std::string& s)
{
  resid_name_ = s;
}
void Residue::set_residue_by_name(const std::string& s)
{
  PhysicalProperty p;
  resid_name_ = s;
  short_name_ = p.get_short_name(s);
  chain_type_ = p.get_chain_type(s);
  mass_ = p.get_mass(s);
  charge_ = p.get_charge(s);
}


inline char Residue::get_chain_ID() const
{
  return chain_ID_;
}
inline void Residue::set_chain_ID(char a)
{
  chain_ID_ = a;
}


inline ChainType Residue::get_chain_type() const
{
  return chain_type_;
}
inline void Residue::set_chain_type(ChainType a)
{
  chain_type_ = a;
}


inline int Residue::get_resid_index() const
{
  return resid_index_;
}
inline void Residue::set_resid_index(int i)
{
  resid_index_ = i;
}


inline int Residue::get_term_flag() const
{
  return term_flag_;
}
inline void Residue::set_term_flag(int i)
{
  // 5: 5';   3: 3';   0: not terminus;
  // -1: N;   1: C;
  term_flag_ = i;
}


inline double Residue::get_resid_charge() const
{
  return charge_;
}
inline void Residue::set_resid_charge(double c)
{
  charge_ = c;
}

inline double Residue::get_resid_mass() const
{
  return mass_;
}
inline void Residue::set_resid_mass(double m)
{
  mass_ = m;
}


inline Atom& Residue::get_atom(int n)
{
  if (atoms_.empty())
  {
    std::cerr << "ERROR: No Atoms in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= int(atoms_.size()))
  {
    std::cerr << "ERROR: Atom index out of range in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return atoms_[n];
}


inline int Residue::add_atom(const Atom& a)
{
  if (a.get_resid_index() != resid_index_)
  {
    return 1;
  }
  for (const Atom& b : atoms_) {
    if (a.get_atom_name() == b.get_atom_name())
      return 0;
  }  // in case of NMR uncertain multi atoms
  atoms_.push_back(a);

  if (a.get_atom_name() == "CA  ")
  {
    C_alpha_ = a;
  }
  if (a.get_atom_name() == "CB  ")
  {
    C_beta_ = a;
  }
  if (a.get_atom_name() == "C3' " || a.get_atom_name() == "S   " || a.get_atom_name() == "DS  ")
  {
    S_ = a;
  }
  if (a.get_atom_name() == "P   " || a.get_atom_name() == "DP  ")
  {
    P_ = a;
  }
  if (a.get_atom_name() == "N1  " || a.get_atom_name() == "B   " || a.get_atom_name() == "DB  ")
  {
    B_ = a;
  }
  if (a.get_atom_flag() == "HETATM" && a.get_element() != "H")
  {
    C_alpha_ = a;
  }
  n_atom_++;
  return 0;
}

inline int Residue::delete_atom(const int i)
{
  if (i >= n_atom_){
    return 1;
  }
  atoms_.erase(atoms_.begin() + i);
  n_atom_--;
  return 0;
}



inline Atom& Residue::get_C_alpha()
{
  if (C_alpha_.get_atom_name() == "") {
    std::cerr << "ERROR: C_alpha not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return C_alpha_;
}

inline Atom& Residue::get_C_beta()
{
  if (C_beta_.get_atom_name() == "") {
    std::cerr << "ERROR: C_beta not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return C_beta_;
}

inline Atom& Residue::get_P()
{
  if (P_.get_atom_name() == "") {
    std::cerr << "ERROR: CG Phosphate not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return P_;
}

inline Atom& Residue::get_S()
{
  if (S_.get_atom_name() == "") {
    std::cerr << "ERROR: CG Sugar not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return S_;
}

inline Atom& Residue::get_B()
{
  if (B_.get_atom_name() == "") {
    std::cerr << "ERROR: CG Base not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return B_;
}

inline void Residue::set_C_alpha()
{
  for (const Atom& b : atoms_) {
    if (b.get_atom_name() == "CA  ")
      C_alpha_ = b;
    break;
  }
}

inline void Residue::set_C_beta()
{
  for (const Atom& b : atoms_) {
    if (b.get_atom_name() == "CB  ")
      C_beta_ = b;
    break;
  }
}

inline void Residue::set_P(const Atom& a)
{
  P_ = a;
}
inline void Residue::set_S(const Atom& a)
{
  S_ = a;
}
inline void Residue::set_B(const Atom& a)
{
  B_ = a;
}

inline void Residue::self_check() const
{
  for (const Atom& a : atoms_) {
    if (a.get_chain_ID() != chain_ID_ || a.get_resid_index() != resid_index_
        || a.get_resid_name() != resid_name_)
    {
      std::cerr << "ERROR: Inconsistent chain ID or residue index or residue type in Residue "
                << resid_index_ << std::endl;
      exit(EXIT_SUCCESS);
    }
  }
}


// void Residue::set_cg_na()
// {
//   int i = 0;
//   Vec3d coor_P(0,0,0);
//   Vec3d coor_C5p(0,0,0),
//       coor_C4p(0,0,0),
//       coor_O4p(0,0,0),
//       coor_C1p(0,0,0),
//       coor_C2p(0,0,0),
//       coor_C3p(0,0,0),
//       coor_O2p(0,0,0),
//       coor_O3p(0,0,0),
//       coor_O5p(0,0,0);
//   Vec3d coor_N1(0,0,0),
//       coor_C2(0,0,0),
//       coor_N3(0,0,0),
//       coor_C4(0,0,0),
//       coor_C5(0,0,0),
//       coor_C6(0,0,0),
//       coor_O2(0,0,0),
//       coor_N4(0,0,0),
//       coor_O4(0,0,0),
//       coor_N6(0,0,0),
//       coor_O6(0,0,0),
//       coor_N2(0,0,0),
//       coor_N7(0,0,0),
//       coor_C8(0,0,0),
//       coor_N9(0,0,0);
//   Vec3d com_P(0,0,0);
//   Vec3d com_S(0,0,0);
//   Vec3d com_B(0,0,0);
//   double mass_C = 12.011;
//   double mass_O = 15.999;
//   double mass_N = 14.001;
//   int n_cs=0;
//   int n_cb=0;
//   int n_os=0;
//   int n_ob=0;
//   int n_nb=0;

//   if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
//   {
//     return;
//   }
//   for (i = 0; i < n_atom_; i++) {
//     std::string aname = atoms_[i].get_atom_name();
//     char c = aname[0];
//     switch (c) {
//       case 'C':
//         if (aname == "C5'") {coor_C5p = atoms_[i].get_coordinates(); n_cs++;}
//         else if (aname == "C1'") {coor_C1p = atoms_[i].get_coordinates(); n_cs++;}
//         else if (aname == "C2'") {coor_C2p = atoms_[i].get_coordinates(); n_cs++;}
//         else if (aname == "C3'") {coor_C3p = atoms_[i].get_coordinates(); n_cs++;}
//         else if (aname == "C4'") {coor_C4p = atoms_[i].get_coordinates(); n_cs++;}
//         else if (aname == "C2 ") { coor_C2 = atoms_[i].get_coordinates(); n_cb++;}
//         else if (aname == "C4 ") { coor_C4 = atoms_[i].get_coordinates(); n_cb++;}
//         else if (aname == "C5 ") { coor_C5 = atoms_[i].get_coordinates(); n_cb++;}
//         else if (aname == "C6 ") { coor_C6 = atoms_[i].get_coordinates(); n_cb++;}
//         else if (aname == "C8 ") { coor_C8 = atoms_[i].get_coordinates(); n_cb++;}
//         break;
//       case 'O':
//         if (aname == "O4'") {coor_O4p = atoms_[i].get_coordinates(); n_os++;}
//         else if (aname == "O2'") {coor_O2p = atoms_[i].get_coordinates(); n_os++;}
//         else if (aname == "O3'") {coor_O3p = atoms_[i].get_coordinates();}
//         else if (aname == "O5'") {coor_O5p = atoms_[i].get_coordinates();}
//         else if (aname == "O2 ") {coor_O2 = atoms_[i].get_coordinates(); n_ob++;}
//         else if (aname == "O4 ") {coor_O4 = atoms_[i].get_coordinates(); n_ob++;}
//         else if (aname == "O6 ") {coor_O6 = atoms_[i].get_coordinates(); n_ob++;}
//         break;
//       case 'N':
//         if (aname == "N1 ") {coor_N1 = atoms_[i].get_coordinates(); n_nb++;}
//         else if (aname == "N2 ") {coor_N2 = atoms_[i].get_coordinates(); n_nb++;}
//         else if (aname == "N3 ") {coor_N3 = atoms_[i].get_coordinates(); n_nb++;}
//         else if (aname == "N4 ") {coor_N4 = atoms_[i].get_coordinates(); n_nb++;}
//         else if (aname == "N6 ") {coor_N6 = atoms_[i].get_coordinates(); n_nb++;}
//         else if (aname == "N7 ") {coor_N7 = atoms_[i].get_coordinates(); n_nb++;}
//         else if (aname == "N9 ") {coor_N9 = atoms_[i].get_coordinates(); n_nb++;}
//         break;
//       default:
//         if (aname == "P  ") coor_P = atoms_[i].get_coordinates();
//     }
//   }
//   com_P = coor_P;
//   com_S = ( (coor_C1p + coor_C2p + coor_C3p + coor_C4p + coor_C5p)
//             * mass_C
//             + ( coor_O2p + coor_O4p ) * mass_O ) * (1/( n_os * mass_O + n_cs * mass_C ));
//   com_B = ( (coor_N1 + coor_N2 + coor_N3 + coor_N4 + coor_N6 + coor_N7 + coor_N9)
//             * mass_N
//             + (coor_O2 + coor_O4 + coor_O6) * mass_O
//             + (coor_C2 + coor_C4 + coor_C5 + coor_C6 + coor_C8)
//             * mass_C )
//           * (1/(n_nb*mass_N + n_ob*mass_O + n_cb*mass_C));
//   P_.set_coords(com_P);
//   if (S_.get_atom_name() != "S  ")
//     S_.set_coords(com_S);
//   if (B_.get_atom_name() != "B  ")
//     B_.set_coords(com_B);
// }


inline int Residue::get_residue_size() const
{
  return n_atom_;
}

// Residue------------------------------------------------------------------
inline Residue::Residue()
{
  resid_name_ = "";
  short_name_ = "0";
  chain_ID_ = -1;
  resid_index_ = -1;
  atoms_.clear();
  n_atom_ = 0;
  charge_ = 0.0;
  mass_ = 100.0;
  term_flag_ = 0;

  C_alpha_.reset();
  C_beta_.reset();
  P_.reset();
  S_.reset();
  B_.reset();
  chain_type_ = none;
}

inline void Residue::reset()
{
  resid_name_ = "";
  short_name_ = "0";
  chain_ID_ = -1;
  resid_index_ = -1;
  atoms_.clear();
  n_atom_ = 0;
  charge_ = 0.0;
  mass_ = 100.0;
  term_flag_ = 0;

  chain_type_ = none;
  P_.reset();
  S_.reset();
  B_.reset();
  C_alpha_.reset();
  C_beta_.reset();
}


// Other functions -----------------------------------------------------------
inline std::ostream& operator<<(std::ostream& o, Residue& r)
{
  int i = 0;
  for (i = 0; i < r.get_residue_size(); i++) {
    o << r.get_atom(i) << std::endl;
  }
  return o;
}

inline double resid_min_distance (Residue& r1, Residue& r2)
{
  int i, j;
  double d = atom_distance(r1.get_atom(0), r2.get_atom(0));  // min_distance;
  double f = 0.0;           // tmp distance;
  for (i = 0; i < r1.get_residue_size(); i++) {
    if (r1.get_atom(i).get_element() == "H")
      continue;
    for (j = 0; j < r2.get_residue_size(); j++) {
      if (r2.get_atom(j).get_element() == "H")
        continue;
      f = atom_distance(r1.get_atom(i), r2.get_atom(j));
      if (d > f)
        d = f;
    }
  }
  return d;
}

inline double resid_ca_distance (Residue& r1, Residue& r2)
{
  double d = -1;           // distance;
  d = atom_distance(r1.get_C_alpha(), r2.get_C_alpha());
  return d;
}

}

#endif

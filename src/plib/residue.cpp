#include "residue.hpp"

namespace pinang {

void Residue::set_resid_name(const std::string& s)
{
  size_t sz = 3;
  resid_name_ = s;
  if (resid_name_.size() < sz)
  {
    resid_name_.resize(sz, ' ');
  }
}
void Residue::set_residue_by_name(const std::string& s)
{
  PhysicalProperty p;
  size_t sz = 3;
  resid_name_ = s;
  if (resid_name_.size() < sz)
  {
    resid_name_.resize(sz, ' ');
  }
  short_name_ = p.get_short_name(s);
  chain_type_ = p.get_chain_type(s);
  mass_ = p.get_mass(s);
  charge_ = p.get_charge(s);
}

Atom& Residue::get_atom(int n)
{
  if (atoms_.empty())
  {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: No Atoms in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= int(atoms_.size()))
  {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: Atom index out of range in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return atoms_[n];
}


int Residue::add_atom(const Atom& a)
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

int Residue::delete_atom(const int i)
{
  if (i >= n_atom_){
    return 1;
  }
  atoms_.erase(atoms_.begin() + i);
  n_atom_--;
  return 0;
}



Atom& Residue::get_C_alpha()
{
  if (C_alpha_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: C_alpha not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return C_alpha_;
}

Atom& Residue::get_C_beta()
{
  if (C_beta_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: C_beta not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return C_beta_;
}

Atom& Residue::get_P()
{
  if (P_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: CG Phosphate not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return P_;
}

Atom& Residue::get_S()
{
  if (S_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: CG Sugar not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return S_;
}

Atom& Residue::get_B()
{
  if (B_.get_atom_name() == "") {
    std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
    std::cerr << "ERROR: CG Base not set in Residue: "
              << resid_index_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return B_;
}

void Residue::set_C_alpha()
{
  for (const Atom& b : atoms_) {
    if (b.get_atom_name() == "CA  ")
      C_alpha_ = b;
    break;
  }
}

void Residue::set_C_beta()
{
  for (const Atom& b : atoms_) {
    if (b.get_atom_name() == "CB  ")
      C_beta_ = b;
    break;
  }
}

void Residue::self_check() const
{
  for (const Atom& a : atoms_) {
    if (a.get_chain_ID() != chain_ID_ || a.get_resid_index() != resid_index_
        || a.get_resid_name() != resid_name_)
    {
      std::cout << " ~               PINANG :: residues.hpp         ~ " << std::endl;
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


// Residue------------------------------------------------------------------
Residue::Residue()
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

void Residue::reset()
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
std::ostream& operator<<(std::ostream& o, Residue& r)
{
  int i = 0;
  int s = r.n_atom_;
  for (i = 0; i < s; i++) {
    o << r.atoms_[i] << std::endl;
  }
  return o;
}

double resid_min_distance(const Residue& r1, const Residue& r2)
{
  int i, j;
  double d = atom_distance(r1.atoms_[0], r2.atoms_[0]);  // min_distance;
  double f = 0.0;           // tmp distance;
  int s1 = r1.n_atom_;
  int s2 = r2.n_atom_;
  for (i = 0; i < s1; i++) {
    if (r1.atoms_[i].get_element() == "H")
      continue;
    for (j = 0; j < s2; j++) {
      if (r2.atoms_[j].get_element() == "H")
        continue;
      f = atom_distance(r1.atoms_[i], r2.atoms_[j]);
      if (d > f)
        d = f;
    }
  }
  return d;
}

double resid_ca_distance(const Residue& r1, const Residue& r2)
{
  double d = -1;           // distance;
  d = atom_distance(r1.C_alpha_, r2.C_alpha_);
  return d;
}

}  // pinang

#include <iomanip>
#include <sstream>
#include "atom.hpp"

namespace pinang {

void Atom::set_atom_name(const std::string& s)
{
  size_t sz = 4;
  atom_name_ = s;
  if (atom_name_.size() < sz)
  {
    atom_name_.resize(sz, ' ');
  }
}

void Atom::set_residue_name(const std::string& s)
{
  size_t sz = 3;
  residue_name_ = s;
  if (residue_name_.size() < sz)
  {
    residue_name_.resize(sz, ' ');
  }
}

// Atom::Atom ==============================================================
Atom::Atom()
{
  record_name_ = "";
  atom_serial_ = 0;
  atom_name_ = "";
  alt_loc_ = ' ';
  residue_name_ = "";
  chain_ID_ = ' ';
  residue_serial_ = 0;
  insert_code_ = ' ';
  coordinate_ = Vec3d(0.0, 0.0, 0.0);
  occupancy_ = 0.0;
  temp_factor_ = 0.0;
  seg_ID_ = "";
  element_ = "";
  charge_ = "";
}

void Atom::reset()
{
  record_name_ = "";
  atom_serial_ = 0;
  atom_name_ = "";
  alt_loc_ = ' ';
  residue_name_ = "";
  chain_ID_ = ' ';
  residue_serial_ = 0;
  insert_code_ = ' ';
  coordinate_ = Vec3d(0.0, 0.0, 0.0);
  occupancy_ = 0.0;
  temp_factor_ = 0.0;
  seg_ID_ = "";
  element_ = "";
  charge_ = "";
}

// outer functions +++++++++++++++++++++++++++++++++++++++++++++
std::ostream& operator<<(std::ostream& o, const Atom& a)
{
  if (a.get_atom_flag() == "ATOM  " || a.get_atom_flag() == "HETATM")
  {
    o << std::setw(6) << a.record_name_
      << std::setw(5) << a.atom_serial_ << " "
      << std::setw(4) << a.atom_name_
      << std::setw(1) << a.alt_loc_
      << std::setw(3) << a.residue_name_ << " "
      << std::setw(1) << a.chain_ID_
      << std::setw(4) << a.residue_serial_
      << std::setw(1) << a.insert_code_ << "   "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
      << std::setw(8) << a.coordinate_.x()
      << std::setw(8) << a.coordinate_.y()
      << std::setw(8) << a.coordinate_.z()
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
      << std::setw(6) << a.occupancy_
      << std::setw(6) << a.temp_factor_ << "      "
      << std::setw(4) << a.seg_ID_
      << std::setw(2) << a.element_
      << std::setw(2) << a.charge_;
  } else {
    o << a.record_name_;
  }
  return o;
}

std::istream& operator>>(std::istream& i, Atom& a)
{
  std::istringstream tmp_sstr;
  std::string pdb_line;

  int tmp_ui;
  std::string tmp_str;
  char tmp_char;
  double tmp_coor_x;
  double tmp_coor_y;
  double tmp_coor_z;
  double tmp_d;

  std::getline(i, pdb_line);
  pdb_line.resize(80, ' ');
  tmp_str = pdb_line.substr(0,6);
  a.record_name_ = tmp_str;

  if (a.get_atom_flag() == "ATOM  " || a.get_atom_flag() == "HETATM")
  {
    tmp_sstr.str(pdb_line.substr(6,5));
    tmp_sstr >> tmp_ui;
    a.atom_serial_ = tmp_ui;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(12,4));
    tmp_sstr >> tmp_str;
    a.set_atom_name(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(16,1));
    tmp_sstr.get(tmp_char);
    a.alt_loc_ = tmp_char;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(17,3));
    tmp_sstr >> tmp_str;
    a.set_residue_name(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(21,1));
    tmp_sstr >> tmp_char;
    a.chain_ID_ = tmp_char;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(22,4));
    tmp_sstr >> tmp_ui;
    a.residue_serial_ = tmp_ui;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(26,1));
    tmp_sstr.get(tmp_char);
    a.insert_code_ = tmp_char;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(30,8));
    tmp_sstr >> tmp_coor_x;
    tmp_sstr.clear();
    tmp_sstr.str(pdb_line.substr(38,8));
    tmp_sstr >> tmp_coor_y;
    tmp_sstr.clear();
    tmp_sstr.str(pdb_line.substr(46,8));
    tmp_sstr >> tmp_coor_z;
    tmp_sstr.clear();
    a.coordinate_ = Vec3d(tmp_coor_x, tmp_coor_y, tmp_coor_z);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(54,6));
    tmp_sstr >> tmp_d;
    a.occupancy_ = tmp_d;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(60,6));
    tmp_sstr >> tmp_d;
    a.temp_factor_ = tmp_d;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(72,4));
    tmp_sstr >> tmp_str;
    a.seg_ID_ = tmp_str;
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(76,2));
    tmp_sstr >> tmp_str;
    a.element_ = tmp_str;
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(78,2));
    tmp_sstr >> tmp_str;
    a.charge_ = tmp_str;
    tmp_sstr.clear();
    tmp_str.clear();
  } else if (a.get_atom_flag() == "MODEL ") {
    tmp_sstr.str(pdb_line.substr(10,4));
    tmp_sstr >> tmp_ui;
    a.atom_serial_ = tmp_ui;  // Actually this is the model index (serial);
    tmp_sstr.clear();
  } else {
    // std::cerr << "ERROR: Wrong PDB ATOM format!" << std::endl;
  }

  if (!i) return i;
  return i;
}

double atom_distance (const Atom& a1, const Atom& a2)
{
  Vec3d v3 = a1.coordinate_ - a2.coordinate_;
  return v3.norm();
}

}  // pinang

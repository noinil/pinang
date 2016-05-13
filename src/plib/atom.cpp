#include "atom.hpp"

namespace pinang {

// std::string Atom::get_atom_flag() const
// {
//   return atom_flag_;
// }
// void Atom::set_atom_flag(const std::string& s)
// {
//   atom_flag_ = s;
// }

// int Atom::get_serial() const
// {
//   return serial_;
// }
// void Atom::set_serial(int i)
// {
//   serial_ = i;
// }

// std::string Atom::get_atom_name() const
// {
//   return atom_name_;
// }
void Atom::set_atom_name(const std::string& s)
{
  size_t sz = 4;
  atom_name_ = s;
  if (atom_name_.size() < sz)
  {
    atom_name_.resize(sz, ' ');
  }
}

// char Atom::get_alt_loc() const
// {
//   return alt_loc_;
// }
// void Atom::set_alt_loc(char a)
// {
//   alt_loc_ = a;
// }

// std::string Atom::get_resid_name() const
// {
//   return resid_name_;
// }
void Atom::set_resid_name(const std::string& s)
{
  size_t sz = 3;
  resid_name_ = s;
  if (resid_name_.size() < sz)
  {
    resid_name_.resize(sz, ' ');
  }
}

// char Atom::get_chain_ID() const
// {
//   return chain_ID_;
// }
// void Atom::set_chain_ID(char a)
// {
//   chain_ID_ = a;
// }

// int Atom::get_resid_index() const
// {
//   return resid_index_;
// }
// void Atom::set_resid_index(int i)
// {
//   resid_index_ = i;
// }

// char Atom::get_icode() const
// {
//   return insert_code_;
// }
// void Atom::set_icode(char a)
// {
//   insert_code_ = a;
// }

// const Vec3d& Atom::get_coordinates() const
// {
//   return coordinate_;
// }
// void Atom::set_coords(const Vec3d& coors)
// {
//   coordinate_ = coors;
// }
// void Atom::set_coords(double x, double y, double z)
// {
//   coordinate_ = Vec3d(x, y, z);
// }

// double Atom::get_occupancy() const
// {
//   return occupancy_;
// }
// void Atom::set_occupancy(double o)
// {
//   occupancy_ = o;
// }

// double Atom::get_temperature_factor() const
// {
//   return temp_factor_;
// }
// void Atom::set_temperature_factor(double f)
// {
//   temp_factor_ = f;
// }

// std::string Atom::get_segment_ID() const
// {
//   return seg_ID_;
// }
// void Atom::set_segment_ID(const std::string& s)
// {
//   seg_ID_ = s;
// }

// std::string Atom::get_element() const
// {
//   return element_;
// }
// void Atom::set_element(const std::string& s)
// {
//   element_ = s;
// }

// std::string Atom::get_charge() const
// {
//   return charge_;
// }
// void Atom::set_charge(const std::string& s)
// {
//   charge_ = s;
// }

// Atom::Atom ==============================================================
Atom::Atom()
{
  atom_flag_ = "";
  serial_ = 0;
  atom_name_ = "";
  alt_loc_ = ' ';
  resid_name_ = "";
  chain_ID_ = ' ';
  resid_index_ = 0;
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
  atom_flag_ = "";
  serial_ = 0;
  atom_name_ = "";
  alt_loc_ = ' ';
  resid_name_ = "";
  chain_ID_ = ' ';
  resid_index_ = 0;
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
    o << std::setw(6) << a.atom_flag_
      << std::setw(5) << a.serial_ << " "
      << std::setw(4) << a.atom_name_
      << std::setw(1) << a.alt_loc_
      << std::setw(3) << a.resid_name_ << " "
      << std::setw(1) << a.chain_ID_
      << std::setw(4) << a.resid_index_
      << std::setw(1) << a.insert_code_ << "   "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
      << std::setw(8) << a.coordinate_.x()
      << std::setw(8) << a.coordinate_.y()
      << std::setw(8) << a.coordinate_.z()
        // << a.coordinate_
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
      << std::setw(6) << a.occupancy_
      << std::setw(6) << a.temp_factor_ << "      "
      << std::setw(4) << a.seg_ID_
      << std::setw(2) << a.element_
      << std::setw(2) << a.charge_;
  } else {
    o << a.atom_flag_;
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
  a.atom_flag_ = tmp_str;

  if (a.get_atom_flag() == "ATOM  " || a.get_atom_flag() == "HETATM")
  {
    tmp_sstr.str(pdb_line.substr(6,5));
    tmp_sstr >> tmp_ui;
    a.serial_ = tmp_ui;
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
    a.set_resid_name(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(21,1));
    tmp_sstr >> tmp_char;
    a.chain_ID_ = tmp_char;
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(22,4));
    tmp_sstr >> tmp_ui;
    a.resid_index_ = tmp_ui;
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
    a.serial_ = tmp_ui;  // Actually this is the model index (serial);
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

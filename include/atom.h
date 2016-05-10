// -*-c++-*-

#ifndef PINANG_ATOM_H_
#define PINANG_ATOM_H_

#include "vec3d.h"

#include <string>
#include <sstream>

namespace pinang {
class Atom
{
 public:
  inline Atom();
  virtual ~Atom() {};

  inline void reset();

  inline std::string get_atom_flag() const;
  inline void set_atom_flag(const std::string& s);

  inline int get_serial() const;
  inline void set_serial(int i);

  inline std::string get_atom_name() const;
  inline void set_atom_name(const std::string& s);

  inline char get_alt_loc() const;
  inline void set_alt_loc(char a);

  inline std::string get_resid_name() const;
  inline void set_resid_name(const std::string& s);

  inline char get_chain_ID() const;
  inline void set_chain_ID(char a);

  inline int get_resid_index() const;
  inline void set_resid_index(int i);

  inline char get_icode() const;
  inline void set_icode(char a);

  inline const Vec3d& get_coordinates() const;
  inline void set_coords(const Vec3d& coors);
  inline void set_coords(double x, double y, double z);

  inline double get_occupancy() const;
  inline void set_occupancy(double o);

  inline double get_temperature_factor() const;
  inline void set_temperature_factor(double f);

  inline std::string get_segment_ID() const;
  inline void set_segment_ID(const std::string& s);

  inline std::string get_element() const;
  inline void set_element(const std::string& s);

  inline std::string get_charge() const;
  inline void set_charge(const std::string& s);

 protected:
  std::string atom_flag_;
  int serial_;
  std::string atom_name_;
  char alt_loc_;
  std::string resid_name_;
  char chain_ID_;
  int resid_index_;
  char insert_code_;
  Vec3d coordinate_;
  double occupancy_;
  double temp_factor_;
  std::string seg_ID_;
  std::string element_;
  std::string charge_;
};

inline std::string Atom::get_atom_flag() const
{
  return atom_flag_;
}
inline void Atom::set_atom_flag(const std::string& s)
{
  atom_flag_ = s;
}

inline int Atom::get_serial() const
{
  return serial_;
}
inline void Atom::set_serial(int i)
{
  serial_ = i;
}

inline std::string Atom::get_atom_name() const
{
  return atom_name_;
}
inline void Atom::set_atom_name(const std::string& s)
{
  size_t sz = 4;
  atom_name_ = s;
  if (atom_name_.size() < sz)
  {
    atom_name_.resize(sz, ' ');
  }
}

inline char Atom::get_alt_loc() const
{
  return alt_loc_;
}
inline void Atom::set_alt_loc(char a)
{
  alt_loc_ = a;
}

inline std::string Atom::get_resid_name() const
{
  return resid_name_;
}
inline void Atom::set_resid_name(const std::string& s)
{
  size_t sz = 3;
  resid_name_ = s;
  if (resid_name_.size() < sz)
  {
    resid_name_.resize(sz, ' ');
  }
}

inline char Atom::get_chain_ID() const
{
  return chain_ID_;
}
inline void Atom::set_chain_ID(char a)
{
  chain_ID_ = a;
}

inline int Atom::get_resid_index() const
{
  return resid_index_;
}
inline void Atom::set_resid_index(int i)
{
  resid_index_ = i;
}

inline char Atom::get_icode() const
{
  return insert_code_;
}
inline void Atom::set_icode(char a)
{
  insert_code_ = a;
}

inline const Vec3d& Atom::get_coordinates() const
{
  return coordinate_;
}
inline void Atom::set_coords(const Vec3d& coors)
{
  coordinate_ = coors;
}
inline void Atom::set_coords(double x, double y, double z)
{
  coordinate_ = Vec3d(x, y, z);
}

inline double Atom::get_occupancy() const
{
  return occupancy_;
}
inline void Atom::set_occupancy(double o)
{
  occupancy_ = o;
}

inline double Atom::get_temperature_factor() const
{
  return temp_factor_;
}
inline void Atom::set_temperature_factor(double f)
{
  temp_factor_ = f;
}

inline std::string Atom::get_segment_ID() const
{
  return seg_ID_;
}
inline void Atom::set_segment_ID(const std::string& s)
{
  seg_ID_ = s;
}

inline std::string Atom::get_element() const
{
  return element_;
}
inline void Atom::set_element(const std::string& s)
{
  element_ = s;
}

inline std::string Atom::get_charge() const
{
  return charge_;
}
inline void Atom::set_charge(const std::string& s)
{
  charge_ = s;
}

// Atom::Atom ==============================================================
inline Atom::Atom()
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

inline void Atom::reset()
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


// other functions +++++++++++++++++++++++++++++++++++++++++++++
inline std::ostream& operator<<(std::ostream& o, const Atom& a)
{
  if (a.get_atom_flag() == "ATOM  " || a.get_atom_flag() == "HETATM")
  {
    o << std::setw(6) << a.get_atom_flag()
      << std::setw(5) << a.get_serial() << " "
      << std::setw(4) << a.get_atom_name()
      << std::setw(1) << a.get_alt_loc()
      << std::setw(3) << a.get_resid_name() << " "
      << std::setw(1) << a.get_chain_ID()
      << std::setw(4) << a.get_resid_index()
      << std::setw(1) << a.get_icode() << "   "
      << a.coordinates()
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
      << std::setw(6) << a.get_occupancy()
      << std::setw(6) << a.get_temperature_factor() << "      "
      << std::setw(4) << a.get_segment_ID()
      << std::setw(2) << a.get_element()
      << std::setw(2) << a.get_charge();
  } else {
    o << a.get_atom_flag();
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
  a.set_atom_flag(tmp_str);

  if (a.get_atom_flag() == "ATOM  " || a.get_atom_flag() == "HETATM")
  {
    tmp_sstr.str(pdb_line.substr(6,5));
    tmp_sstr >> tmp_ui;
    a.set_serial(tmp_ui);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(12,4));
    tmp_sstr >> tmp_str;
    a.set_atom_name(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(16,1));
    tmp_sstr.get(tmp_char);
    a.set_alt_loc(tmp_char);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(17,3));
    tmp_sstr >> tmp_str;
    a.set_resid_name(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(21,1));
    tmp_sstr >> tmp_char;
    a.set_chain_ID(tmp_char);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(22,4));
    tmp_sstr >> tmp_ui;
    a.set_resid_index(tmp_ui);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(26,1));
    tmp_sstr.get(tmp_char);
    a.set_icode(tmp_char);
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
    a.set_coords(tmp_coor_x, tmp_coor_y, tmp_coor_z);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(54,6));
    tmp_sstr >> tmp_d;
    a.set_occupancy(tmp_d);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(60,6));
    tmp_sstr >> tmp_d;
    a.set_temperature_factor(tmp_d);
    tmp_sstr.clear();

    tmp_sstr.str(pdb_line.substr(72,4));
    tmp_sstr >> tmp_str;
    a.set_segment_ID(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(76,2));
    tmp_sstr >> tmp_str;
    a.set_element(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();

    tmp_sstr.str(pdb_line.substr(78,2));
    tmp_sstr >> tmp_str;
    a.set_charge(tmp_str);
    tmp_sstr.clear();
    tmp_str.clear();
  } else if (a.get_atom_flag() == "MODEL ") {
    tmp_sstr.str(pdb_line.substr(10,4));
    tmp_sstr >> tmp_ui;
    a.set_serial(tmp_ui);  // Actually this is the model index (serial);
    tmp_sstr.clear();
  } else {
    // std::cerr << "ERROR: Wrong PDB ATOM format!" << std::endl;
  }

  if (!i) return i;
  return i;
}

inline double atom_distance (Atom& a1, Atom& a2)
{
  Vec3d v3 = a1.coordinates() - a2.coordinates();
  return v3.norm();
}

}
#endif

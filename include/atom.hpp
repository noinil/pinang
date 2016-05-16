/*!
************************************************************
@file atom.hpp
@brief Definition of class ATOM.

In this file class ATOM is defined.  The data format is defined according to the
PDB style.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 17:53
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_ATOM_H_
#define PINANG_ATOM_H_

#include "vec3d.hpp"

#include <string>

namespace pinang {

class Atom
{
 public:
  Atom();
  virtual ~Atom() {};

  void reset();

  std::string get_atom_flag() const { return atom_flag_; }
  void set_atom_flag(const std::string& s) { atom_flag_ = s; }

  int get_serial() const { return serial_; }
  void set_serial(int i) { serial_ = i; }

  std::string get_atom_name() const { return atom_name_; }
  void set_atom_name(const std::string&);

  char get_alt_loc() const { return alt_loc_; }
  void set_alt_loc(char a) { alt_loc_ = a; }

  std::string get_resid_name() const { return resid_name_; }
  void set_resid_name(const std::string&);

  char get_chain_ID() const { return chain_ID_; }
  void set_chain_ID(char a) { chain_ID_ = a; }

  int get_resid_index() const { return resid_index_; }
  void set_resid_index(int i) { resid_index_ = i; }

  char get_icode() const { return insert_code_; }
  void set_icode(char a) { insert_code_ = a; }

  const Vec3d& get_coordinates() const { return coordinate_; }
  void set_coords(const Vec3d& coors) { coordinate_ = coors; }
  void set_coords(double x, double y, double z) { coordinate_ = Vec3d(x, y, z); }

  double get_occupancy() const { return occupancy_; }
  void set_occupancy(double o) { occupancy_ = o; }

  double get_temperature_factor() const { return temp_factor_; }
  void set_temperature_factor(double f) { temp_factor_ = f; }

  std::string get_segment_ID() const { return seg_ID_; }
  void set_segment_ID(const std::string& s) { seg_ID_ = s; }

  std::string get_element() const { return element_; }
  void set_element(const std::string& s) { element_ = s; }

  std::string get_charge() const { return charge_; }
  void set_charge(const std::string& s) { charge_ = s; }

  friend std::ostream& operator<<(std::ostream&, const Atom&);
  friend std::istream& operator>>(std::istream&, Atom&);
  friend double atom_distance (const Atom&, const Atom&);
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

double atom_distance (const Atom&, const Atom&);
}
#endif

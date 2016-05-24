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

/*!
************************************************************
@brief Physical atomistic type, with properties in the PDB format.

Read in PDB lines and store the information together within this class.
************************************************************
*/
class Atom
{
 public:
  // ************************************************************
  //! @brief Create an "empty" Atom object.
  //! @return An Atom object.
  // ************************************************************
  Atom();
  virtual ~Atom() {};

  // ************************************************************
  //! @brief Reset properties of Atom.
  // ************************************************************
  void reset();

  // ************************************************************
  //! @brief Get the "ATOM  " record keyword in PDB.
  //! @return Record word.
  // ************************************************************
  std::string get_atom_flag() const { return record_name_; }
  // ************************************************************
  //! @brief Set the "ATOM  " record keyword of atom.
  //! @param Record word.
  // ************************************************************
  void set_atom_flag(const std::string& s) { record_name_ = s; }

  // ************************************************************
  //! @brief Get serial number of atom.
  //! @return Serial number of atom.
  // ************************************************************
  int get_serial() const { return serial_; }
  // ************************************************************
  //! @brief Set serial number of atom.
  //! @param Serial number of atom.
  // ************************************************************
  void set_serial(int i) { serial_ = i; }

  // ************************************************************
  //! @brief Get atom name.
  //! @return Atom name.
  // ************************************************************
  std::string get_atom_name() const { return atom_name_; }
  // ************************************************************
  //! @brief Set atom name.
  //! @param Atom name.
  // ************************************************************
  void set_atom_name(const std::string&);

  // ************************************************************
  //! @brief Get alternate location indicator.
  //! @return Alternate location indicator.
  // ************************************************************
  char get_alt_loc() const { return alt_loc_; }
  // ************************************************************
  //! @brief Set alternate location indicator.
  //! @param Alternate location indicator.
  // ************************************************************
  void set_alt_loc(char a) { alt_loc_ = a; }

  // ************************************************************
  //! @brief Get residue name.
  //! @return Residue name.
  // ************************************************************
  std::string get_resid_name() const { return resid_name_; }
  // ************************************************************
  //! @brief Set residue name.
  //! @param Residue name.
  // ************************************************************
  void set_resid_name(const std::string&);

  // ************************************************************
  //! @brief Get chain identifier.
  //! @return Chain identifier.
  // ************************************************************
  char get_chain_ID() const { return chain_ID_; }
  // ************************************************************
  //! @brief Set chain identifier.
  //! @param Chain identifier, such as 'A', 'B', 'X'...
  // ************************************************************
  void set_chain_ID(char a) { chain_ID_ = a; }

  // ************************************************************
  //! @brief Get residue serial number.
  //! @return Residue serial number.
  // ************************************************************
  int get_resid_index() const { return resid_index_; }
  // ************************************************************
  //! @brief Set residue serial number.
  //! @param Residue serial number.
  // ************************************************************
  void set_resid_index(int i) { resid_index_ = i; }

  // ************************************************************
  //! @brief Get Code for insertion of residues.
  //! @return Code for insertion of residues.
  // ************************************************************
  char get_icode() const { return insert_code_; }
  // ************************************************************
  //! @brief Set Code for insertion of residues.
  //! @param Code for insertion of residues.
  // ************************************************************
  void set_icode(char a) { insert_code_ = a; }

  // ************************************************************
  //! @brief Get coordinates of atom.
  //! @return Coordinates of atom.
  // ************************************************************
  const Vec3d& get_coordinates() const { return coordinate_; }
  // ************************************************************
  //! @brief Get coordinates of atom from Vec3d type coordinates.
  //! @param Vec3d type coordinates.
  // ************************************************************
  void set_coords(const Vec3d& coors) { coordinate_ = coors; }
  // ************************************************************
  //! @brief Get coordinates of atom from three real numbers.
  //! @param x, y, and z, three orthogonal coordinates.
  // ************************************************************
  void set_coords(double x, double y, double z) { coordinate_ = Vec3d(x, y, z); }

  // ************************************************************
  //! @brief Get occupancy.
  //! @return Occupancy.
  // ************************************************************
  double get_occupancy() const { return occupancy_; }
  // ************************************************************
  //! @brief Set occupancy.
  //! @param Occupancy.
  // ************************************************************
  void set_occupancy(double o) { occupancy_ = o; }

  // ************************************************************
  //! @brief Get temperature factor.
  //! @return Temperature factor.
  // ************************************************************
  double get_temperature_factor() const { return temp_factor_; }
  // ************************************************************
  //! @brief Set temperature factor.
  //! @param Temperature factor.
  // ************************************************************
  void set_temperature_factor(double f) { temp_factor_ = f; }

  // ************************************************************
  //! @brief Get segment identifier.
  //! @return Segment identifier.
  // ************************************************************
  std::string get_segment_ID() const { return seg_ID_; }
  // ************************************************************
  //! @brief Set segment identifier.
  //! @param Segment identifier.
  // ************************************************************
  void set_segment_ID(const std::string& s) { seg_ID_ = s; }

  // ************************************************************
  //! @brief Get element symbol.
  //! @return Element symbol.
  // ************************************************************
  std::string get_element() const { return element_; }
  // ************************************************************
  //! @brief Set element symbol.
  //! @param Element symbol.
  // ************************************************************
  void set_element(const std::string& s) { element_ = s; }

  // ************************************************************
  //! @brief Get atom charge.
  //! @return Atom charge.
  // ************************************************************
  std::string get_charge() const { return charge_; }
  // ************************************************************
  //! @brief Set atom charge.
  //! @param Atom charge.
  // ************************************************************
  void set_charge(const std::string& s) { charge_ = s; }

  // ************************************************************
  //! @brief Output PDB format information of Atom.
  // ************************************************************
  friend std::ostream& operator<<(std::ostream&, const Atom&);
  // ************************************************************
  //! @brief Read in PDB information to Atom.
  // ************************************************************
  friend std::istream& operator>>(std::istream&, Atom&);
  friend double atom_distance (const Atom&, const Atom&);
 protected:
  std::string record_name_;  //!< Atom flag from PDB.
  int serial_;             //!< Atom serial number in PDB.
  std::string atom_name_;  //!< Atom name in PDB.
  char alt_loc_;           //!< Alternate location indicator from PDB.
  std::string resid_name_; //!< Residue name of the atom in PDB.
  char chain_ID_;          //!< Chain identifier in PDB.
  int resid_index_;        //!< Residue sequence number in PDB.
  char insert_code_;       //!< Code for insertion of residues from PDB.
  Vec3d coordinate_;       //!< Orthogonal coordinates (x, y, z).
  double occupancy_;       //!< Occupancy from PDB.
  double temp_factor_;     //!< Temperature factor in PDB.
  std::string seg_ID_;     //!< Segment identifier in PDB.
  std::string element_;    //!< Element symbol in PDB.
  std::string charge_;     //!< Charge on the atom.
};

// ************************************************************
//! @brief Compute distance between two Atoms.
//! @param Two Atom objects.
//! @return Distance.
// ************************************************************
double atom_distance (const Atom&, const Atom&);
}
#endif

/*!
************************************************************
@file particle.hpp
@brief Definition of class Particle.

In this file class Particle is defined.  Particle can be an atom, a residue, or
a cluster of many different basic molecule types, depending on the model.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:06
@copyright GNU Public License V3.0
************************************************************
*/


#ifndef PINANG_PARTICLE_H_
#define PINANG_PARTICLE_H_

#include "vec3d.hpp"

#include <string>

namespace pinang {

/*!
************************************************************
@brief Whatever basic mass point.

This class is designed to represent any mass point in molecular dynamics
simulations.  The important thing are the mass and charge, which is used in the
topology file.
************************************************************
*/
class Particle
{
 public:
  // ************************************************************
  //! @brief Create an "empty" Particle object.
  //! @return A Particle object.
  // ************************************************************
  Particle();
  virtual ~Particle() {};

  // ************************************************************
  //! @brief Reset properties of Particle.
  // ************************************************************
  void reset();

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
  //! @brief Get particle charge.
  //! @return Particle charge.
  // ************************************************************
  double get_charge() const { return charge_; }
  // ************************************************************
  //! @brief Set particle charge.
  //! @param Particle charge.
  // ************************************************************
  void set_charge(double c) { charge_ = c; }

  // ************************************************************
  //! @brief Get particle mass.
  //! @return Particle mass.
  // ************************************************************
  double get_mass() const { return mass_; }
  // ************************************************************
  //! @brief Set particle mass.
  //! @param Particle mass.
  // ************************************************************
  void set_mass(double m) { mass_ = m; }

  // ************************************************************
  //! @brief Read in information to Particle.
  // ************************************************************
  friend std::istream& operator>>(std::istream&, Particle&);

 protected:
  std::string atom_name_;   //!< Atom name of particle.
  std::string resid_name_;  //!< Residue name of particle.
  int resid_index_;         //!< Residue sequence number of particle.
  double charge_;           //!< Charge of particle.
  double mass_;             //!< Mass of particle.
};

}
#endif

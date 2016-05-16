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
  Particle();
  virtual ~Particle() {};

  void reset();

  std::string get_atom_name() const { return atom_name_; }
  void set_atom_name(const std::string&);

  std::string get_resid_name() const { return resid_name_; }
  void set_resid_name(const std::string&);

  int get_resid_index() const { return resid_index_; }
  void set_resid_index(int i) { resid_index_ = i; }

  double get_charge() const { return charge_; }
  void set_charge(double c) { charge_ = c; }

  double get_mass() const { return mass_; }
  void set_mass(double m) { mass_ = m; }

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

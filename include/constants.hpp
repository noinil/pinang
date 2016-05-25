/*!
  @file constants.hpp
  @brief Definition of basic constants.

  In this file a list of basic constants such as unit_of_mass, spring_constant for
  energy functions, cutoff, chain_type, and residue_name are defined.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 17:58
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_CONSTANTS_H_
#define PINANG_CONSTANTS_H_

#include <map>
#include <string>
#include <iostream>

namespace pinang {

extern double g_cutoff;  //!< Cutoff for atomistic distances.

const long double k_pi = 3.14159265358979323846;  //!< @f$\pi@f$.

// Units.
const double k_u_mass = 1.0;  //!< Unit of mass.

// CG MD Parameters.
const double k_K_bond = 100.0;      //!< Energy function parameter for bonds.
const double k_K_angle = 20.0;      //!< Energy function parameter for angles.
const double k_K_dihedral_1 = 1.0;  //!< Energy function parameter for dihedral angles.
const double k_K_dihedral_3 = 0.5;  //!< Energy function parameter for dihedral angles (optional).
const double k_K_native = 1.0;      //!< Energy function parameter for native contacts.
const double k_K_nonnative = 1.0;   //!< Energy function parameter for non-native contacts.

//! Chain chemical types.
enum ChainType {none=0, protein=1, DNA=2, RNA=3, water=4, ion=5, other=6, na=7};

/*!
  @brief Some of the physical properties of biomolecules.

  A collection of maps that can translate residue name into charges, masses and
  chemical chain types.
*/
class PhysicalProperty
{
 public:
  //! @brief Create an "empty" PhysicalProperty object.
  PhysicalProperty();
  ~PhysicalProperty();

  //! @brief Translate residue name into short name.
  //! @param Residue name.
  //! @return Short name.
  std::string get_short_name(const std::string&);
  //! @brief Translate residue name into charge.
  //! @param Residue name.
  //! @return Charge.
  double get_charge(const std::string&);
  //! @brief Translate residue name into mass.
  //! @param Residue name.
  //! @return Mass.
  double get_mass(const std::string&);
  //! @brief Translate residue name into chain type.
  //! @param Residue name.
  //! @return Chain type.
  ChainType get_chain_type(const std::string&);

 private:
  std::map<std::string, std::string> map_resName_shortName;  //!< Mapping of residue name to short name.
  std::map<std::string, double> map_resName_charge;          //!< Mapping of residue name to residue charge.
  std::map<std::string, double> map_resName_mass;            //!< Mapping of residue name to residue mass.
  std::map<std::string, ChainType> map_resName_chainType;    //!< Mapping of residue name to chain type.
};

}

#endif

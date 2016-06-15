/*!
  @file residue.hpp
  @brief Definition of class Residue.

  In this file class Residue is defined.  Residue contains a list of atoms which
  are connected in a line.  Other than the basic physical atoms, residues have
  psudo atoms such as C-alpha beads.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 18:09
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_RESIDUE_H_
#define PINANG_RESIDUE_H_

#include <vector>

#include "atom.hpp"
#include "constants.hpp"

namespace pinang {

/*!
  @brief Biomolecular building block, made of a few atoms.

  A few atoms comprises a residue.  The class Residue records the coordinates of
  atoms, the atom names, residue name, residue serial number and other properties.

  @todo initialization or settings of CG particles in Residue (cg_P_, cg_S_, cg_B_).
*/
class Residue
{
 public:
  //! @brief Create an "empty" Residue object.
  //! @return An Residue object.
  Residue();
  virtual ~Residue() {v_atoms_.clear();}

  //! @brief Reset properties of Residue.
  void reset();

  //! @brief Get residue name.
  //! @return Residue name.
  std::string get_residue_name() const { return residue_name_; }
  //! @brief Get short residue name.
  //! @return Short residue name.
  std::string get_short_name() const { return short_name_; }
  //! @brief Set residue name.
  //! @param Residue name.
  void set_residue_name(const std::string&);
  //! @brief Set residue physical properties by residue name.
  //! @param Residue name.
  void set_residue_by_name(const std::string&);

  //! @brief Get chain identifier.
  //! @return Chain identifier.
  char get_chain_ID() const { return chain_ID_; }
  //! @brief Set chain identifier.
  //! @param Chain identifier, such as 'A', 'B', 'X'...
  void set_chain_ID(char a) { chain_ID_ = a; }

  //! @brief Get chain type.
  //! @return Chain type. (enum type)
  ChainType get_chain_type() const { return chain_type_; }
  //! @brief Set chain type.
  //! @param Chain type.
  void set_chain_type(ChainType a) { chain_type_ = a; }

  //! @brief Get residue serial number.
  //! @return Residue serial number.
  int get_residue_serial() const { return residue_serial_; }
  //! @brief Set residue serial number.
  //! @param Residue serial number.
  void set_residue_serial(int i) { residue_serial_ = i; }

  //! @brief Get terminus information.
  //! @return Terminus information.
  //! @retval 5: 5' of DNA (RNA).
  //! @retval 3: 3' of DNA (RNA).
  //! @retval -1: N-terminus of protein.
  //! @retval 1: C-terminus of protein.
  //! @retval 0: Not terminus.
  int get_terminus_flag() const { return terminus_flag_; }
  //! @brief Set terminus information.
  //! @param Terminus flag.
  void set_terminus_flag(int i) { terminus_flag_ = i; }

  //! @brief Get Residue charge.
  //! @return Residue charge.
  double get_residue_charge() const { return charge_; }
  //! @brief Set Residue charge.
  //! @param Residue charge.
  void set_residue_charge(double c) { charge_ = c; }

  //! @brief Get Residue mass.
  //! @return Residue mass.
  double get_residue_mass() const { return mass_; }
  //! @brief Set Residue mass.
  //! @param Residue mass.
  void set_residue_mass(double m) { mass_ = m; }

  //! @brief Self check before additional operations.
  void self_check() const;

  //! @brief Get a Atom object from Residue.
  //! @param Serial number of the Atom.
  //! @return Atom.
  Atom& get_atom(int);
  //! @brief Add an Atom object to Residue.
  //! @param Atom.
  //! @return Status of adding atom to residue.
  //! @retval 0: Success.
  int add_atom(const Atom&);
  //! @brief Deleting an Atom object from Residue.
  //! @param Serial number of Atom.
  //! @return Status of deleting atom from residue.
  //! @retval 0: Success.
  //! @retval 1: Failure.
  int delete_atom(const int);

  //! @brief Get residue size.
  //! @return Residue size.
  int get_size() const { return n_atom_; }

  //! @brief Get CG particle @f$C_\alpha@f$.
  //! @return CG particle @f$C_\alpha@f$.
  Atom& get_cg_C_alpha();
  //! @brief Get CG particle @f$C_\beta@f$.
  //! @return CG particle @f$C_\beta@f$.
  Atom& get_cg_C_beta();
  //! @brief Get CG particle P (phosphate).
  //! @return CG particle P (phosphate).
  Atom& get_cg_P();
  //! @brief Get CG particle S (sugar).
  //! @return CG particle S (sugar).
  Atom& get_cg_S();
  //! @brief Get CG particle B (base).
  //! @return CG particle B (base).
  Atom& get_cg_B();
  //! @brief Set CG particle @f$C_\alpha@f$.
  void set_cg_C_alpha();
  //! @brief Set CG particle @f$C_\beta@f$.
  void set_cg_C_beta();
  //! @brief Set CG particle P.
  //! @param CG Atom phosphate.
  void set_cg_P(const Atom& a) { cg_P_ = a; }
  //! @brief Set CG particle S.
  //! @param CG Atom sugar.
  void set_cg_S(const Atom& a) { cg_S_ = a; }
  //! @brief Set CG particle B.
  //! @param CG Atom base.
  void set_cg_B(const Atom& a) { cg_B_ = a; }

  //! @brief Output PDB format information of Atom.
  friend std::ostream& operator<<(std::ostream&, Residue&);
  friend double residue_min_distance(const Residue&, const Residue&);
  friend double residue_min_distance(const Residue&, const Residue&, Atom&, Atom&);
  friend double residue_ca_distance(const Residue&, const Residue&);
 protected:
  std::string residue_name_;   //!< Residue name from PDB.
  std::string short_name_;   //!< Short name of residue.
  char chain_ID_;            //!< Chain identifier from PDB.
  int residue_serial_;          //!< Residue sequence number in PDB.
  std::vector<Atom> v_atoms_;  //!< A set of atom objects.
  int n_atom_;               //!< Number of atoms in Residue.
  double charge_;            //!< Charge of Residue.
  double mass_;              //!< Mass of Residue.

  Atom cg_C_alpha_;             //!< CG particle @f$C_\alpha@f$.
  Atom cg_C_beta_;              //!< CG particle @f$C_\beta@f$.
  Atom cg_P_;                   //!< CG particle P (phosphate).
  Atom cg_S_;                   //!< CG particle S (sugar).
  Atom cg_B_;                   //!< CG particle B (base).
  ChainType chain_type_;     //!< Chain type.

  // 5: 5'; 3: 3'; -1: N;  1: C;  0: not terminus;
  int terminus_flag_;            //!< Flag indicating whether the Residue is at terminus.
};

//! @brief Compute the minimal distance between two Residues.
//! @param Two Residue objects.
//! @return Distance.
double residue_min_distance(const Residue&, const Residue&);
//! @brief Compute the minimal distance between two Residues.
//! @param Two Residue objects, two atom objects to store min-distance paires.
//! @return Distance.
double residue_min_distance(const Residue&, const Residue&, Atom&, Atom&);
//! @brief Compute the @f$C_\alpha@f$ distance between two Residues.
//! @param Two Residue objects.
//! @return @f$C_\alpha@f$ Distance.
double residue_ca_distance(const Residue&, const Residue&);
}

#endif

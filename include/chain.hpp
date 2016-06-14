/*!
  @file chain.hpp
  @brief Definition of class Chain.

  In this file class Chain is defined.  Chain is biologically a collection of
  tandemly connected residues.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 17:56
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_CHAIN_H_
#define PINANG_CHAIN_H_

#include "residue.hpp"

namespace pinang {

/*!
  @brief Biological molecule chains consisting of residues.

  This class represents biomolecular chains, composing a chain of residues, either
  amino acids or nucleic acids.   Read and store information from PDB files.

  @todo initialization or settings of CG particles in DNA residues (cg_P_, cg_S_, cg_B_).
*/
class Chain
{
 public:
  //! @brief Create an "empty" Chain object.
  //! @return An Chain object.
  Chain();
  virtual ~Chain() {v_residues_.clear();}

  //! @brief Reset properties of Chain.
  void reset();

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

  //! @brief Get a Residue object from Chain.
  //! @param Serial number of the Residue.
  //! @return Residue.
  Residue& get_residue(int);
  //! @brief Add a Residue object to Chain.
  //! @param Residue.
  //! @return Status of adding residues to chain.
  //! @retval 0: Success.
  int add_residue(const Residue&);

  //! @brief Get chain length (number of residues included).
  //! @return Chain length.
  int get_size() const { return n_residue_; }

  //! @brief Print sequence of the Chain.
  //! @param Option to output short style (1) or full name (3).
  void output_sequence(int) const;
  //! @brief Output sequence information to a fasta-style file.
  void output_sequence_fasta(std::ostream&, std::string) const;
  //! @brief Self check before additional operations.
  void self_check();

  //! @brief Output PDB of CG beads.
  void output_cg_pdb(std::ostream&, int&, int&);
  //! @brief Output coordinates of CG beads.
  void output_cg_crd(std::ostream&);

  //! @brief Output physical properties to topology file.
  void output_top_mass(std::ostream&, int&, int&);
  //! @brief Output bonded interactions to topology file.
  void output_top_bond(std::ostream&, int&, int&);
  //! @brief Output angle interactions to topology file.
  void output_top_angle(std::ostream&, int&, int&);
  //! @brief Output dihedral angle interactions to topology file.
  void output_top_dihedral(std::ostream&, int&, int&);

  //! @brief Output bonded interactions to forcefield parm file.
  void output_ffparm_bond(std::ostream&, int&);
  //! @brief Output angle interactions to forcefield parm file.
  void output_ffparm_angle(std::ostream&, int&);
  //! @brief Output dihedral angle interactions to forcefield parm file.
  void output_ffparm_dihedral(std::ostream&, int&);
  //! @brief Output non-bonded native contact interactions to forcefield parm file.
  void output_ffparm_native(std::ostream&);

  //! @brief Get native contact number intra-chain.
  //! @return Native contact number.
  int get_native_contact_number();

  //! @brief Add two chains together (connect two chains).
  //! @param Two chains.
  //! @return A new chain.
  friend Chain operator+(const Chain&, const Chain&);

  //! @brief Output PDB format information of Chain.
  friend std::ostream& operator<<(std::ostream&, Chain&);

 protected:
  char chain_ID_;                  //!< Chain identifier in PDB.
  ChainType chain_type_;           //!< Chain chemical composition.
  int n_residue_;                  //!< Number of residues.
  std::vector<Residue> v_residues_;  //!< A collection of residue objects.
};

}

#endif

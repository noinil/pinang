/*!
  @file model.hpp
  @brief Definition of class Model.

  In this file class Model is defined.  Model is a configuration of a biomolecule,
  usually including one or more chains.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 18:05
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_MODEL_H_
#define PINANG_MODEL_H_

#include "chain.hpp"

namespace pinang {

/*!
  @brief Certain configuration of molecules stored in PDB.

  A model is usually a certain conformation of biomolecules stored in PDB
  structures (especially in NMR structures).  A model consists of one or several
  chains.
*/
class Model
{
 public:
  //! @brief Create an "empty" Model object.
  //! @return A Model object.
  Model();
  virtual ~Model() {v_chains_.clear();}

  //! @brief Reset properties of Model.
  void reset();

  //! @brief Get model identifier.
  //! @return Model identifier.
  int get_model_serial() const { return model_serial_; }
  //! @brief Set model identifier.
  //! @param Model identifier.
  void set_model_serial(int n) { model_serial_ = n; }

  //! @brief Get a Chain object from Model.
  //! @param Serial number of the Chain.
  //! @return Chain.
  Chain& get_chain(unsigned int);
  //! @brief Add a Chain object to Model.
  //! @param Chain.
  //! @return Status of adding Chain to Model.
  //! @retval 0: Success.
  void add_chain(Chain&);

  //! @brief Print sequence of the Model.
  //! @param Option to output short style (1) or full name (3).
  void output_sequence(int) const;
  //! @brief Output sequence information to a fasta-style file.
  void output_sequence_fasta(std::ostream&, std::string) const;

  //! @brief Get number of chains in the Model.
  //! @return Number of chains in the Model.
  int get_size() const { return n_chain_; }

  //! @brief Output PDB of CG beads.
  void output_cg_pdb(std::ostream&);
  //! @brief Output coordinates of CG beads.
  void output_cg_crd(std::ostream&);

  //! @brief Output physical properties to topology file.
  void output_top_mass(std::ostream&);
  //! @brief Output bonded interactions to topology file.
  void output_top_bond(std::ostream&);
  //! @brief Output angle interactions to topology file.
  void output_top_angle(std::ostream&);
  //! @brief Output dihedral angle interactions to topology file.
  void output_top_dihedral(std::ostream&);

  //! @brief Output bonded interactions to forcefield parm file.
  void output_ffparm_bond(std::ostream&);
  //! @brief Output angle interactions to forcefield parm file.
  void output_ffparm_angle(std::ostream&);
  //! @brief Output dihedral angle interactions to forcefield parm file.
  void output_ffparm_dihedral(std::ostream&);
  //! @brief Output non-bonded native contact interactions to forcefield parm file.
  void output_ffparm_nonbonded(std::ostream&);

  //! @brief Output statistics of protein-DNA pairwise residue-residue distances.
  void output_statistics_pro_DNA_contact_pairs(std::ostream&);

  //! @brief Output PDB format information of Chain.
  friend std::ostream& operator<<(std::ostream&, Model&);

 protected:
  int model_serial_;               //!< Model serial number.
  std::vector<Chain> v_chains_;  //!< A set of chains objects.
  int n_chain_;                //!< Number of chains in Model.
};

}

#endif

/*!
************************************************************
@file PDB.hpp
@brief Definition of class PDB.

This file defines class PDB, which contains structure to store information
from PDB structure files and provides functions to print sequence, re-print
clear PDB coordinates.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 17:00
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_PDB_H
#define PINANG_PDB_H

#include "model.hpp"

namespace pinang{

/*!
************************************************************
@brief Read in biomolecular structures from pdb file.

Store the molecules (including coordinates, residue names, etc.) and output
sequence or coordinates.
************************************************************
*/

class PDB
{
 public:
  // ************************************************************
  //! @brief Create a PDB object by read from PDB flie.
  //! @param PDB file name.
  //! @return A PDB object.
  // ************************************************************
  PDB(const std::string& s);
  virtual ~PDB() {v_models_.clear();}

  // ************************************************************
  //! @brief Get PDB file name.
  //! @return PDB file name.
  // ************************************************************
  std::string get_pdb_name() const { return PDB_file_name_; }

  // ************************************************************
  //! @brief Get access to a model in a PDB structure.
  //! @param Index of the model.
  //! @return A Model object.
  // ************************************************************
  Model& get_model(unsigned int);

  // ************************************************************
  //! @brief Get the number of models in PDB.
  //! @return Number of models.
  // ************************************************************
  int get_size() const { return n_model_; }

  // ************************************************************
  //! @brief Print sequence of whole PDB.
  //! @param Option to output short style (1) or full name (3).
  // ************************************************************
  void output_sequence(int) const;

  // ************************************************************
  //! @brief Output sequence information to a fasta-style file.
  // ************************************************************
  void output_sequence_fasta(std::ostream&) const;

  // ************************************************************
  //! @brief Output PDB format information to ostream.
  // ************************************************************
  friend std::ostream& operator<<(std::ostream&, PDB&);

 protected:
  std::string PDB_file_name_;  //!< PDB flie name.
  std::vector<Model> v_models_;  //!< A collection of model objects in PDB file.
  int n_model_;                //!< Number of models in PDB flie.
};

}

#endif

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

#include <iostream>
#include <string>

namespace pinang{

class PDB
{
 public:
  PDB(const std::string& s);
  virtual ~PDB() {models_.clear();}

  std::string get_pdb_name() const { return PDB_file_name_; }

  Model& get_model(unsigned int);
  int get_n_models() const { return n_model_; }

  void print_sequence(int) const;
  void output_fasta(std::ostream&) const;

  friend std::ostream& operator<<(std::ostream&, PDB&);

 protected:
  std::string PDB_file_name_;  //!< PDB flie name.
  std::vector<Model> models_;  //!< A collection of model objects in PDB file.
  int n_model_;                //!< Number of models in PDB flie.
};

}

#endif

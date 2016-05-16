/*!
************************************************************
@file model.hpp
@brief Definition of class Model.

In this file class Model is defined.  Model is a configuration of a biomolecule,
usually including one or more chains.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:05
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_MODEL_H_
#define PINANG_MODEL_H_

#include "chain.hpp"

namespace pinang {

/*!
************************************************************
@brief Certain configuration of molecules stored in PDB.

A model is usually a certain conformation of biomolecules stored in PDB
structures (especially in NMR structures).  A model consists of one or several
chains.
************************************************************
*/
class Model
{
 public:
  Model();
  virtual ~Model() {chains_.clear();}

  void reset();

  int get_model_ID() const { return model_ID_; }
  void set_model_ID(int n) { model_ID_ = n; }

  Chain& get_chain(unsigned int);
  void add_chain(Chain&);

  void print_sequence(int) const;
  void output_fasta(std::ostream&, std::string) const;

  int get_model_size() const { return n_chain_; }

  void output_cg_pos(std::ostream&);
  void output_top_mass(std::ostream&);
  void output_top_bond(std::ostream&);
  void output_top_angle(std::ostream&);
  void output_top_dihedral(std::ostream&);
  void output_top_nonbonded(std::ostream&);

  friend std::ostream& operator<<(std::ostream&, Model&);

 protected:
  int model_ID_;               //!< Model serial number.
  std::vector<Chain> chains_;  //!< A set of chains objects.
  int n_chain_;                //!< Number of chains in Model.
};

}

#endif

// -*-c++-*-

#ifndef PINANG_PDB_H
#define PINANG_PDB_H

#include "model.h"

#include <fstream>
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
  std::string PDB_file_name_;
  std::vector<Model> models_;
  int n_model_;
};

}

#endif

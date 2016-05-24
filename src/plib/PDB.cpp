#include <fstream>
#include "PDB.hpp"

namespace pinang {

PDB::PDB(const std::string& s)
{
  Atom atom_tmp;
  Residue resid_tmp;
  Chain chain_tmp;
  Model model_tmp;

  PDB_file_name_ = s;
  n_model_ = 0;
  v_models_.clear();

  std::ifstream ifile(PDB_file_name_.c_str());
  if (!ifile.is_open())
  {
    std::cout << " ~             PINANG :: PDB.h                ~ " << std::endl;
    std::cerr << " ERROR: Cannot read file: " << s << std::endl;
    exit(EXIT_FAILURE);
  }

  while (ifile.good()) {
    ifile >> atom_tmp;
    if (ifile.fail())
    {
      break;
    }

    if (atom_tmp.get_atom_flag() == "MODEL ")
    {
      model_tmp.reset();
      model_tmp.set_model_serial(atom_tmp.get_serial());

      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "TER   ")
    {
      if (resid_tmp.get_residue_size() != 0)
      {
        chain_tmp.add_residue(resid_tmp);
        chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
        chain_tmp.set_chain_type(resid_tmp.get_chain_type());
      }
      model_tmp.add_chain(chain_tmp);

      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "ENDMDL")
    {
      if (resid_tmp.get_residue_size() != 0)
      {
        chain_tmp.add_residue(resid_tmp);
        chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
        chain_tmp.set_chain_type(resid_tmp.get_chain_type());
      }
      if (chain_tmp.get_chain_length() != 0)
      {
        model_tmp.add_chain(chain_tmp);
      }
      v_models_.push_back(model_tmp);
      n_model_++;

      model_tmp.reset();
      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "END   ")
    {
      if (resid_tmp.get_residue_size() != 0)
      {
        chain_tmp.add_residue(resid_tmp);
        chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
        chain_tmp.set_chain_type(resid_tmp.get_chain_type());
      }
      if (chain_tmp.get_chain_length() != 0)
      {
        model_tmp.add_chain(chain_tmp);
      }
      if (model_tmp.get_model_size() != 0)
      {
        v_models_.push_back(model_tmp);
        n_model_++;
      }

      model_tmp.reset();
      chain_tmp.reset();
      resid_tmp.reset();
      atom_tmp.reset();
    }
    if (atom_tmp.get_atom_flag() == "ATOM  " )
    {
      if (resid_tmp.add_atom(atom_tmp))
      {
        if (resid_tmp.get_residue_size() != 0)
        {
          chain_tmp.add_residue(resid_tmp);
          if (resid_tmp.get_atom(0).get_atom_flag() == "HETATM")
          {
            chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
            chain_tmp.set_chain_type(resid_tmp.get_chain_type());
            model_tmp.add_chain(chain_tmp);
            chain_tmp.reset();
          }

          resid_tmp.reset();
        }
        resid_tmp.set_residue_by_name(atom_tmp.get_resid_name());
        resid_tmp.set_chain_ID(atom_tmp.get_chain_ID());
        resid_tmp.set_residue_serial(atom_tmp.get_residue_serial());

        resid_tmp.add_atom(atom_tmp);
      }
    }
    if (atom_tmp.get_atom_flag() == "HETATM")
    {
      if (resid_tmp.add_atom(atom_tmp))
      {
        if (resid_tmp.get_residue_size() != 0)
        {
          chain_tmp.add_residue(resid_tmp);
          chain_tmp.set_chain_ID(resid_tmp.get_chain_ID());
          chain_tmp.set_chain_type(resid_tmp.get_chain_type());
          model_tmp.add_chain(chain_tmp);

          chain_tmp.reset();
          resid_tmp.reset();
        }
        resid_tmp.set_residue_serial(atom_tmp.get_residue_serial());
        resid_tmp.set_resid_name(atom_tmp.get_resid_name());
        resid_tmp.set_chain_ID(atom_tmp.get_chain_ID());
        resid_tmp.add_atom(atom_tmp);
      }
    }
  }

  ifile.close();
}

Model& PDB::get_model(unsigned int n)
{
  if (v_models_.empty())
  {
    std::cout << " ~             PINANG :: PDB.h                ~ " << std::endl;
    std::cerr << "ERROR: No Model found in this PDB: "
              << PDB_file_name_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  if (n >= v_models_.size())
  {
    std::cout << " ~             PINANG :: PDB.h                ~ " << std::endl;
    std::cerr << "ERROR: Model number out of range in PDB: "
              << PDB_file_name_ << std::endl;
    exit(EXIT_SUCCESS);
  }
  return v_models_[n];
}

void PDB::print_sequence(int n) const
{
  if (n != 1 && n != 3)
  {
    std::cerr << " Usage: PINANG::PDB.print_sequence(): \n"
              << "       n = 1: 1-char residue name;\n"
              << "       n = 3: 3-char residue name.\n"
              << std::endl;
    exit(EXIT_SUCCESS);
  }
  v_models_[0].print_sequence(n);
}

void PDB::output_fasta(std::ostream & f_fasta) const
{
  std::string s = PDB_file_name_;
  for (int i = 0; i < 4; i++) {
    s.pop_back();
  }

  v_models_[0].output_fasta(f_fasta, s);
}

std::ostream& operator<<(std::ostream& o, PDB& p)
{
  int i = 0;
  int s = p.n_model_;
  for (i = 0; i < s; i++) {
    o << p.v_models_[i] << std::endl;
  }
  return o;
}

}

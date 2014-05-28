// -*-c++-*-

#ifndef PINANG_PDB_H
#define PINANG_PDB_H

#include "model.h"

#include <iostream>
#include <fstream>
#include <string>

namespace pinang{

    class PDB
    {
    public:
        PDB(const std::string& s);
        virtual ~PDB() {_models.clear();}

        inline std::string pdb_name() const;

        inline int output();



    protected:
        std::string _PDB_file_name;
        std::vector<Model> _models;
        int _n_model;
    };

    PDB::PDB(const std::string& s)
    {
        Atom atom_tmp;
        Residue resid_tmp;
        Chain chain_tmp;
        Model model_tmp;
        unsigned int tmp_ri = 0;
        _PDB_file_name = s;

        std::ifstream ifile(_PDB_file_name.c_str());
        if (!ifile.is_open())
        {
            std::cerr << "ERROR: Cannot read file: " << s << std::endl;
            std::cerr << "Program terminating." << std::endl;
            exit(EXIT_FAILURE);
        }

        // std::cout << "Reading from PDB file: " << s << std::endl;

        while (ifile.good()) {
            ifile >> atom_tmp;
            // std::cout << atom_tmp << std::endl;
            if (ifile.fail())
            {
                break;
            }

            if (atom_tmp.atom_flag() == "MODEL ")
            {
                model_tmp.reset();
                chain_tmp.reset();
                resid_tmp.reset();
                atom_tmp.reset();
                tmp_ri = 0;
            }
            if (atom_tmp.atom_flag() == "TER   ")
            {
                chain_tmp.add_residue(resid_tmp);
                chain_tmp.set_chain_ID(resid_tmp.chain_ID());
                model_tmp.add_chain(chain_tmp);
                chain_tmp.reset();
                resid_tmp.reset();
                tmp_ri = 0;
            }
            if (atom_tmp.atom_flag() == "ENDMDL")
            {
                _models.push_back(model_tmp);
                _n_model++;
                model_tmp.reset();
                chain_tmp.reset();
                resid_tmp.reset();
                atom_tmp.reset();
                tmp_ri = 0;
            }
            if (atom_tmp.atom_flag() == "END   ")
            {
                // std::cout << "model_tmp.size(): " << model_tmp.m_model_size() << std::endl;
                if (model_tmp.m_model_size() != 0)
                {
                    _models.push_back(model_tmp);
                    _n_model++;
                }
                model_tmp.reset();
                chain_tmp.reset();
                resid_tmp.reset();
                atom_tmp.reset();
                tmp_ri = 0;
            }
            if (atom_tmp.atom_flag() == "ATOM  ")
            {
                if (resid_tmp.add_atom(atom_tmp))
                {
                    if (tmp_ri != 0)
                    {
                        chain_tmp.add_residue(resid_tmp);
                        resid_tmp.reset();
                    }
                    tmp_ri = atom_tmp.resid_index();
                    resid_tmp.set_resid_name(atom_tmp.resid_name());
                    resid_tmp.set_chain_ID(atom_tmp.chain_ID());
                    resid_tmp.set_resid_index(atom_tmp.resid_index());

                    resid_tmp.add_atom(atom_tmp);
                }
            }
        }

        ifile.close();
    }

    inline std::string PDB::pdb_name() const
    {
        return _PDB_file_name;
    }

    inline int PDB::output()
    {
        std::cout << "PDB entry: " << _PDB_file_name << std::endl;
        int i = 0;
        for (i = 0; i < _models.size(); i++) {
            std::cout << "MODEL: " << i << std::endl;
            std::cout << _models[i] << std::endl;
        }
        return 0;
    }

}

#endif

// -*-c++-*-

#ifndef PINANG_PDB_H
#define PINANG_PDB_H

#include "model.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

namespace pinang{

    class PDB
    {
    public:
        PDB(const std::string& s);
        virtual ~PDB() {_models.clear();}

        inline std::string pdb_name() const;

        inline Model& m_model(unsigned int n);
        inline int n_models() const;

        int native_contact_number() const;

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
        _n_model = 0;
        _models.clear();

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

    inline Model& PDB::m_model(unsigned int n)
    {
        if (_models.empty())
        {
            std::cerr << "ERROR: No Model found in this PDB: "
                      << _PDB_file_name << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _models.size())
            {
                std::cerr << "ERROR: Model number out of range in PDB: "
                          << _PDB_file_name << std::endl;
                exit(EXIT_SUCCESS);
            } else {
                return _models[n];
            }
        }
    }

    inline std::string PDB::pdb_name() const
    {
        return _PDB_file_name;
    }
    inline int PDB::n_models() const
    {
        return _n_model;
    }

    inline std::ostream& operator<<(std::ostream& o, PDB& p)
    {
        int i = 0;
        for (i = 0; i < p.n_models(); i++) {
            o << "MODEL: " << i << std::endl;
            o << p.m_model(i) << std::endl;
        }
        return o;
    }


    int PDB::native_contact_number() const
    {
        int nc = 0;             // native contact number;
        int i = 0;
        for (i = 0; i < _n_model ; i++) {

        }
        return nc;
    }

}

#endif

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

        inline void print_sequence(int n) const;

        // void contact_map(double c);

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
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~               PINANG :: PDB                ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << " ERROR: Cannot read file: " << s << std::endl;
            std::cerr << " Program terminating." << std::endl;
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
                tmp_ri = atom_tmp.serial();
                model_tmp.set_model_ID(tmp_ri);

                tmp_ri = 0;
                chain_tmp.reset();
                resid_tmp.reset();
                atom_tmp.reset();
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
                if (resid_tmp.m_residue_size() != 0)
                {
                    chain_tmp.add_residue(resid_tmp);
                    chain_tmp.set_chain_ID(resid_tmp.chain_ID());
                    model_tmp.add_chain(chain_tmp);
                }

                _models.push_back(model_tmp); // push back tmp model;
                _n_model++;

                model_tmp.reset();
                chain_tmp.reset();
                resid_tmp.reset();
                atom_tmp.reset();
                tmp_ri = 0;
            }
            if (atom_tmp.atom_flag() == "END   ")
            {
                if (resid_tmp.m_residue_size() != 0)
                {
                    chain_tmp.add_residue(resid_tmp);
                    chain_tmp.set_chain_ID(resid_tmp.chain_ID());
                    model_tmp.add_chain(chain_tmp);
                }

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
            if (atom_tmp.atom_flag() == "ATOM  " || atom_tmp.atom_flag() == "HETATM")
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
                    if (atom_tmp.atom_flag() == "HETATM")
                    {
                        resid_tmp.set_chain_ID('+');
                    } else {
                        resid_tmp.set_chain_ID(atom_tmp.chain_ID());
                    }
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
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~               PINANG :: PDB                ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << "ERROR: No Model found in this PDB: "
                      << _PDB_file_name << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _models.size())
            {
                std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
                std::cout << " ~               PINANG :: PDB                ~ " << std::endl;
                std::cout << " ============================================== " << std::endl;
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

    /*             _       _
    //  _ __  _ __(_)_ __ | |_   ___  ___  __ _
    // | '_ \| '__| | '_ \| __| / __|/ _ \/ _` |
    // | |_) | |  | | | | | |_  \__ \  __/ (_| |
    // | .__/|_|  |_|_| |_|\__| |___/\___|\__, |
    // |_|                                   |_|
    */
    inline void PDB::print_sequence(int n) const
    {
        if (n != 1 && n != 3)
        {
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~               PINANG :: PDB                ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << " Usage: PINANG::PDB.print_sequence(): \n"
                 << "       n = 1: 1-char residue name;\n"
                 << "       n = 3: 3-char residue name.\n"
                 << std::endl;
            exit(EXIT_SUCCESS);
        }
        _models[0].print_sequence(n);
        // std::cout << std::endl;
    }


    inline std::ostream& operator<<(std::ostream& o, PDB& p)
    {
        int i = 0;
        for (i = 0; i < p.n_models(); i++) {
            o << p.m_model(i) << std::endl;
        }
        return o;
    }

}

#endif

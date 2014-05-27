// -*-c++-*-

#ifndef PINANG_PDB_READ_H
#define PINANG_PDB_READ_H

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
        unsigned int tmp_ri = 0;

        ifstream ifile(s.c_str());
        if (!ifile.is_open())
        {
            std::cerr << "ERROR: Cannot read file: " << infilename << std::endl;
            std::cerr << "Program terminating." << endl;
            exit(EXIT_FAILURE);
        }
        while (ifile.good()) {
            ifile >> atom_tmp;
        }
    }
}

#endif

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "model.h"

using namespace std;

int main(int argc, char *argv[])
{
    pinang::Model model;
    pinang::Chain chain_tmp;
    pinang::Residue resid_tmp;
    pinang::Atom atom_tmp;

    std::string infilename = argv[1];
    ifstream ifile(infilename.c_str());
    std::cout << "Reading from PDB " << infilename << std::endl;
    if (!ifile.is_open())
    {
        std::cerr << "Could not open file: " << infilename << std::endl;
        std::cerr << "Program terminating." << endl;
        exit(EXIT_FAILURE);
    }

    std::string pdb_line;
    std::string atom_flag;
    std::istringstream tmp_sstr;

    while (ifile.good()) {
        // getline( ifile, pdb_line);
        ifile >> atom_tmp;
        if (ifile.fail())
        {
            break;
        }
        // pdb_line.resize(80, ' ');
        // atom_flag = pdb_line.substr(0,6);

        // if (atom_flag == "ATOM  ")
        // {
        //     tmp_sstr.str (pdb_line);
        //     tmp_sstr >> atom_tmp;

        //     std::cout << atom_tmp << std::endl;
        //     tmp_sstr.clear();
        // }
        if (atom_tmp.atom_name() != "NULL")
        {
            std::cout << atom_tmp << std::endl;
        }
    }

    return 0;
}

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
    std::stringstream tmp_sstr;
    unsigned int _tmp_serial;
    std::string _tmp_atom_name;
    std::string _tmp_resid_name;

    while (ifile.good()) {
        getline( ifile, pdb_line);
        if (ifile.fail())
        {
            break;
        }
        pdb_line.resize(80, ' ');
        atom_flag = pdb_line.substr(0,6);
        // std::cout << atom_flag << std::endl;
        if (atom_flag == "ATOM  ")
        {
            tmp_sstr.str ( pdb_line.substr(6,5));
            tmp_sstr >> _tmp_serial;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(12,4));
            tmp_sstr >> _tmp_atom_name;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(17,3));
            tmp_sstr >> _tmp_resid_name;
            tmp_sstr.clear();

            std::cout << "serial number: " << setw(4) << _tmp_serial
                      << " atom name: " << setw(4) << _tmp_atom_name
                      << " resid name: " << setw(4) << _tmp_resid_name
                      << std::endl;
        }
    }

    return 0;
}

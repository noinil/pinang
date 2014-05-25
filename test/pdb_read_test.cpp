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

    unsigned int _tmp_serial;
    std::string _tmp_atom_name;
    char _tmp_alt_loc;
    std::string _tmp_resid_name;
    char _tmp_chain_ID;
    unsigned int _tmp_resid_index;
    char _tmp_insert_code;
    pinang::Vec3d _tmp_coordinates;
    double _tmp_occupancy;
    double _tmp_temp_factor;
    std::string _tmp_seg_ID;
    std::string _tmp_element;
    std::string _tmp_charge;

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

            tmp_sstr.str ( pdb_line.substr(16,1));
            tmp_sstr.get(_tmp_alt_loc);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(17,3));
            tmp_sstr >> _tmp_resid_name;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(21,1));
            tmp_sstr >> _tmp_chain_ID;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(22,4));
            tmp_sstr >> _tmp_resid_index;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(26,1));
            tmp_sstr.get(_tmp_insert_code);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(30,24));
            tmp_sstr >> _tmp_coordinates;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(54,6));
            tmp_sstr >> _tmp_occupancy;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(60,6));
            tmp_sstr >> _tmp_temp_factor;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(72,4));
            tmp_sstr >> _tmp_seg_ID;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(76,2));
            tmp_sstr >> _tmp_element;
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(78,2));
            tmp_sstr >> _tmp_charge;
            tmp_sstr.clear();

            atom_tmp.set_serial(_tmp_serial);
            atom_tmp.set_atom_name(_tmp_atom_name);
            atom_tmp.set_alt_loc(_tmp_alt_loc);
            atom_tmp.set_resid_name(_tmp_resid_name);
            atom_tmp.set_chain_ID(_tmp_chain_ID);
            atom_tmp.set_resid_index(_tmp_resid_index);
            atom_tmp.set_icode(_tmp_insert_code);
            atom_tmp.set_coords(_tmp_coordinates);
            atom_tmp.set_occupancy(_tmp_occupancy);
            atom_tmp.set_temperature_factor(_tmp_temp_factor);
            atom_tmp.set_segment_ID(_tmp_seg_ID);
            atom_tmp.set_element(_tmp_element);
            atom_tmp.set_charge(_tmp_charge);

            std::cout << atom_tmp << std::endl;

            // std::cout << "ATOM  " << setw(5) << _tmp_serial << " "
            //           << setw(4) << _tmp_atom_name << " "
            //           << setw(3) << _tmp_resid_name << " "
            //           << setw(1) << _tmp_chain_ID
            //           << setw(4) << _tmp_resid_index << "    "
            //           << _tmp_coordinates
            //           << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
            //           << setw(6) << _tmp_occupancy
            //           << setw(6) << _tmp_temp_factor << "      "
            //           << setw(4) << _tmp_seg_ID
            //           << setw(2) << _tmp_element
            //           << setw(2) << _tmp_charge
            //           << std::endl;
        }
    }

    return 0;
}

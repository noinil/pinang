#include <iostream>
#include <fstream>

#include "PDB.h"

using namespace std;

int main(int argc, char *argv[])
{
    std::string infilename = argv[1];
    std::string outfilename = argv[2];

    pinang::PDB pdb1(infilename);
    std::ofstream ofile(outfilename.c_str());

    // ------------------ PDB output test;
    // std::cout << pdb1 << std::endl;

    // -------------------- PDB seq test;
    // pdb1.print_sequence(1);
    // pdb1.print_sequence(3);
    // pdb1.print_sequence(2);

    // ------------------ PDB output cg pos
    pdb1.m_model(0).output_ca_pos(ofile);

    return 0;
}

#include <iostream>
#include <fstream>

#include "PDB.h"

using namespace std;

int main(int argc, char *argv[])
{
    std::string infilename = argv[1];

    pinang::PDB pdb1(infilename);
    // std::cout << pdb1 << std::endl;

    pdb1.print_sequence(1);
    pdb1.print_sequence(3);
    pdb1.print_sequence(2);

    // pdb1.contact_map();

    return 0;
}

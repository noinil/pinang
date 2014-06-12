#include "PDB.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
    std::string infilename = argv[1];

    pinang::PDB pdb1(infilename);

    std::cout << " ============================================== "
              << std::endl
              << " Sequence of PDB "
              << pdb1.pdb_name()
              << " :"
              << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cout << " 1-char-aa-name : ----------------------------- " << std::endl;
    pdb1.print_sequence(1);
    std::cout << " ---------------------------------------------- " << std::endl;
    std::cout << " 3-char-aa-name : ----------------------------- " << std::endl;
    pdb1.print_sequence(3);
    std::cout << " ---------------------------------------------- " << std::endl;

    return 0;
}

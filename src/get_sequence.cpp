#include "PDB.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~           PINANG sequence print            ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;

    if (argc < 2 || argc >=3)
    {
        std::cerr << " ERROR: too few or too many arguments!" << std::endl;
        std::cerr << " Usage: "
                  << argv[0]
                  << " some.pdb" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string infilename = argv[1];

    pinang::PDB pdb1(infilename);

    std::cout << " Sequence of PDB "
              << pdb1.pdb_name()
              << " :"
              << std::endl;
    std::cout << std::endl;
    std::cout << " 1-char-aa-name : ----------------------------- " << std::endl;
    pdb1.print_sequence(1);
    std::cout << " ---------------------------------------------- " << std::endl;
    std::cout << " 3-char-aa-name : ----------------------------- " << std::endl;
    pdb1.print_sequence(3);
    std::cout << " ---------------------------------------------- " << std::endl;

    return 0;
}

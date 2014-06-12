#include "PDB.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~         PINANG PDB output/extract          ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    // std::cerr << " ERROR: too few or too many arguments!" << std::endl;
    std::cout << " Usage: p_pdb_extract some.pdb" << std::endl;
    // exit(EXIT_FAILURE);

    std::string infilename = argv[1];

    pinang::PDB pdb1(infilename);

    int i = pdb1.n_models();

    if (i > 1)
    {
        std::cout << " Please input Module No. :" << std::endl;
        std::cin >> i;
        std::cout << pdb1.m_model(i) << std::endl;
    } else {
        std::cout << pdb1 << std::endl;
    }

    return 0;
}

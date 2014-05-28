#include <iostream>
#include <fstream>

#include "PDB.h"

using namespace std;

int main(int argc, char *argv[])
{
    std::string infilename = argv[1];

    pinang::PDB pdb1(infilename);
    pdb1.output();

    return 0;
}

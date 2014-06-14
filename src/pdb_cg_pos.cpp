#include "PDB.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~      PINANG PDB CG positions output        ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cout << " Usage: "
              << argv[0]
              << " -f some.pdb [-o output.pos] [-m module]" << std::endl;
    std::cout << " ============================================== " << std::endl;

    int opt, mod_index = 0;
    int mod_flag = 0;
    int in_flag = 0;

    std::string infilename = "some.pdb";
    std::string outfilename = "out.pos";

    while ((opt = getopt(argc, argv, "o:m:f:")) != -1) {
        switch (opt) {
        case 'o':
            outfilename = optarg;
            break;
        case 'm':
            mod_index = atoi(optarg);
            mod_flag = 1;
            break;
        case 'f':
            infilename = optarg;
            in_flag = 1;
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o output.pos] [-m module]" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o output.pos] [-m module]" << std::endl;
        exit(EXIT_FAILURE);
    }
    pinang::PDB pdb1(infilename);

    std::ofstream ofile(outfilename.c_str());
    if (mod_flag == 1) {
        std::cout << " Extracting C-alpha coordinates of MODULE " << mod_index
                  << " of " << infilename
                  << " to " << outfilename
                  << std::endl;
        pdb1.m_model(mod_index - 1).output_ca_pos(ofile);
    } else {
        std::cout << " Extracting C-alpha coordinates of MODULE "
                  << pdb1.m_model(0).model_ID()
                  << " of " << infilename
                  << " to " << outfilename
                  << std::endl;
        pdb1.m_model(0).output_ca_pos(ofile);
    }

    return 0;
}

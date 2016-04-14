#include "PDB.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~               PINANG PDB output/extract                ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt, mod_index = 0;
    int mod_flag = 0;
    int in_flag = 0;

    std::string infilename = "some.pdb";
    std::string outfilename = "out.pdb";

    while ((opt = getopt(argc, argv, "o:m:f:h")) != -1) {
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
        case 'h':
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o output.pdb] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o output.pdb] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o output.pdb] [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    std::ofstream ofile(outfilename.c_str());
    if (mod_flag == 1) {
        std::cout << " Extracting MODULE " << mod_index
                  << " of " << infilename
                  << " to " << outfilename
                  << std::endl;
        ofile << pdb1.m_model(mod_index - 1) << std::endl;
    } else if (pdb1.n_models() > 1) {
        std::cout << " Please choose a MODULE: " ;
        std::cin >> mod_index;
        std::cout << " Extracting MODULE " << mod_index
                  << " of " << infilename
                  << " to " << outfilename
                  << std::endl;
        ofile << pdb1.m_model(mod_index - 1) << std::endl;
    } else {
        std::cout << " Re-output " << infilename
                  << " to " << outfilename
                  << std::endl;
        ofile << pdb1 << std::endl;
    }

    return 0;
}

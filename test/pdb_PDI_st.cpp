#include "PDB.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>

using namespace std;

int main(int argc, char *argv[])
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // whatever

    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~                  PINANG DNA curvature                  ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt, mod_index = 0;
    int mod_flag = 0;
    int in_flag = 0;

    std::string infilename = "some.pdb";

    std::string axis_name = "_axis.pdb";
    std::string back_name = "_backbone.pdb";
    std::string out_name = "_curve.dat";

    while ((opt = getopt(argc, argv, "o:x:b:m:f:h")) != -1) {
        switch (opt) {
        case 'o':
            out_name = optarg;
            break;
        case 'x':
            axis_name = optarg;
            break;
        case 'b':
            back_name = optarg;
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
                      << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                      << " [-b _backbone.pdb] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                      << " [-b _backbone.pdb] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " ERROR: need parameter for option -f: " << std::endl
                  << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                  << " [-b _backbone.pdb] [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    // std::ofstream back_file(back_name.c_str());
    // std::ofstream axis_file(axis_name.c_str());
    // std::ofstream out_file(out_name.c_str());

    if (mod_flag != 1) {
        if (pdb1.n_models() == 1)
        {
            mod_index = 1;
        } else {
            std::cout << " Please choose a MODULE: " ;
            std::cin >> mod_index;
        }
    }

    std::cout << " Analyzing DNA curvature of MODULE " << mod_index
              << " of " << infilename  << " ... "
              << std::endl;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // main
    //  _ _                 _      _ _
    // ( | )_ __ ___   __ _(_)_ __( | )
    //  V V| '_ ` _ \ / _` | | '_ \V V
    //     | | | | | | (_| | | | | |
    //     |_| |_| |_|\__,_|_|_| |_|
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    pinang::Chain c0;
    for (int i = 0; i < pdb1.m_model(mod_index - 1).m_model_size(); i++) {
        c0 = c0 + pdb1.m_model(mod_index - 1).m_chain(i);
    }
    for (int i = 0; i < c0.m_chain_length(); i++) {
        std::cout << i << "  "
                  << c0.m_residue(i).chain_type() << "  "
                  << c0.m_residue(i).resid_name()
                  << std::endl;
    }


    return 0;
}

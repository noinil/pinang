#include "PDB.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~           PINANG sequence print            ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;

    std::string infilename = "some.pdb";
    std::string out_name = "seq.fasta";

    int opt;
    int in_flag = 0;

    while ((opt = getopt(argc, argv, "o:f:h")) != -1) {
        switch (opt) {
        case 'o':
            out_name = optarg;
            break;
        case 'f':
            infilename = optarg;
            in_flag = 1;
            break;
        case 'h':
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o seq.fasta] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o seq.fasta] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " ERROR: need parameter for option -f: " << std::endl
                  << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o seq.fasta] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    std::ofstream out_file(out_name.c_str());

    std::cout << " Sequence of PDB "
              << pdb1.pdb_name()
              << " :"
              << std::endl;
    std::cout << " Total number of chains: "
              << pdb1.m_model(0).m_model_size()
              << std::endl;
    std::cout << std::endl;
    std::cout << " 1-char-aa-name : ----------------------------- " << std::endl;
    pdb1.print_sequence(1);
    std::cout << " ---------------------------------------------- " << std::endl;
    std::cout << " 3-char-aa-name : ----------------------------- " << std::endl;
    pdb1.print_sequence(3);
    std::cout << " ---------------------------------------------- " << std::endl;
    pdb1.output_fasta(out_file);

    return 0;
}

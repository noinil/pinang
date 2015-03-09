#include "PDB.h"

#include <iostream>
#include <vector>
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

    std::vector<int> pro_res_id;
    std::vector<int> dna_res_id;

    std::string infilename = "some.pdb";
    std::string outfilename = "chains.info";

    while ((opt = getopt(argc, argv, "m:f:h")) != -1) {
        switch (opt) {
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
    std::ofstream pro_file("pro.pdb");
    std::ofstream dna_file("dna.pdb");

    if (mod_flag == 1) {
        pinang::Model mdl = pdb1.m_model(mod_index);
    } else if (pdb1.n_models() > 1) {
        std::cout << " Please choose a MODULE: " ;
        std::cin >> mod_index;
        std::cout << " Extracting MODULE " << mod_index
                  << " of " << infilename
                  << " to " << outfilename
                  << std::endl;
        pinang::Model mdl = pdb1.m_model(mod_index);
    } else {
        std::cout << " Extracting MODULE " << 0
                  << " of " << infilename
                  << std::endl;
        pinang::Model mdl = pdb1.m_model(0);
        for (int i = 0; i < mdl.m_model_size(); i++) {
            if (mdl.m_chain(i).chain_type() == 1) {
                pro_res_id.push_back(i);
            } else if (mdl.m_chain(i).chain_type() == 2 or mdl.m_chain(i).chain_type() == 7) {
                dna_res_id.push_back(i);
            } else {
                continue;
            }
        }

        // -------------------- output pro.pdb --------------------
        for (unsigned int i = 0; i < pro_res_id.size(); i++) {
            int cid = pro_res_id[i];
            pinang::Chain chn = mdl.m_chain(cid);
            pro_file << chn;
        }
        pro_file << "END" << std::endl;

        // -------------------- output dna.pdb --------------------
        for (unsigned int i = 0; i < dna_res_id.size(); i++) {
            int cid = dna_res_id[i];
            pinang::Chain chn = mdl.m_chain(cid);
            chn.set_chain_type(pinang::DNA);
            // std::cout << i << " " << chn.chain_type() << std::endl;
            dna_file << chn;
        }
        dna_file << "END" << std::endl;

        int cc1b, cc1e, cc2b, cc2e;
        cc1b = 1;
        cc1e = cc1b + dna_res_id.size() - 1;
        cc2b = 1 + cc1e;
        cc2e = cc2b + pro_res_id.size() - 1;
        ofile << "<<<< unit_and_state" << std::endl;
        ofile << "i_seq_read_style = 1" << std::endl;
        ofile << "i_go_native_read_style = 1" << std::endl;
        if (cc1b < cc1e) {
            ofile << cc1b << "-" << cc1e << "   dna2   dna.pdb" << std::endl;
            std::cout << cc1b << "-" << cc1e << "   dna2   dna.pdb" << std::endl;
        }
        else {
            ofile << cc1b << "   dna2   dna.pdb" << std::endl;
            std::cout << cc1b << "   dna2   dna.pdb" << std::endl;
        }
        if (cc2b < cc2e) {
            ofile << cc2b << "-" << cc2e << "   protein   pro.pdb" << std::endl;
            std::cout << cc2b << "-" << cc2e << "   protein   pro.pdb" << std::endl;
        }
        else {
            ofile << cc2b << "   protein   pro.pdb" << std::endl;
            std::cout << cc2b << "   protein   pro.pdb" << std::endl;
        }
        ofile << ">>>>\n" << std::endl;

        ofile << "<<<< energy_function" << std::endl;
        if (cc1b < cc1e)
            ofile << "LOCAL(" << cc1b
                  << "-" << cc1e << ")    L_DNA2" << std::endl;
        else
            ofile << "LOCAL(" << cc1b
                  << ")    L_DNA2" << std::endl;
        if (cc2b < cc2e)
            ofile << "LOCAL(" << cc2b
                  << "-" << cc2e << ")    L_AICG2_PLUS" << std::endl;
        else
            ofile << "LOCAL(" << cc2b
                  << ")    L_AICG2_PLUS" << std::endl;
        if (cc1b < cc1e)
            ofile << "NLOCAL(" << cc1b
                  << "-" << cc1e << "/"
                  << cc1b << "-" << cc1e
                  << ")    DNA2 ELE" << std::endl;
        else
            ofile << "NLOCAL(" << cc1b << "/"
                  << cc1b << ")    DNA2 ELE" << std::endl;
        if (cc2b < cc2e)
            ofile << "NLOCAL(" << cc2b
                  << "-" << cc2e << "/"
                  << cc2b << "-" << cc2e
                  << ")    AICG2 EXV ELE" << std::endl;
        else
            ofile << "NLOCAL(" << cc2b << "/"
                  << cc2b << ")    AICG2 EXV ELE" << std::endl;
        ofile <<  "NLOCAL(" << cc1b << "-" << cc1e << "/"
              << cc2b << "-" << cc2e << ")    EXV ELE" << std::endl;
    }

    return 0;
}

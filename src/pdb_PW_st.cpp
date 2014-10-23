#include "PDB.h"
#include "vec3d.h"
#include "constants.h"

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

    std::string out_name = "PW_dist.dat";

    while ((opt = getopt(argc, argv, "o:m:f:h")) != -1) {
        switch (opt) {
        case 'o':
            out_name = optarg;
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
                      << " -f some.pdb [-o PW_dist.dat] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o PW_dist.dat] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " ERROR: need parameter for option -f: " << std::endl
                  << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o PW_dist.dat] [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    std::ofstream out_file(out_name.c_str());

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
    pinang::Residue r0, r1;
    pinang::chain_t cti1;
    for (int i = 0; i < pdb1.m_model(mod_index - 1).m_model_size(); i++) {
        if (pdb1.m_model(mod_index - 1).m_chain(i).chain_type() != pinang::water)
            continue;
        for (int m = 0; m < pdb1.m_model(mod_index - 1).m_chain(i).m_chain_length(); m++) {
            r0 = pdb1.m_model(mod_index - 1).m_chain(i).m_residue(m);

            double min_dist_W_PRO = 0; // min distance water -- protein
            double min_dist_W_A = 0; // min distance water -- DNA
            double min_dist_W_T = 0; // min distance water -- DNA
            double min_dist_W_G = 0; // min distance water -- DNA
            double min_dist_W_C = 0; // min distance water -- DNA
            double min_dist_W_S = 0; // min distance water -- DNA
            double min_dist_W_P = 0; // min distance water -- DNA

            for (int j = 0; j < pdb1.m_model(mod_index - 1).m_model_size(); j++) {
                if (pdb1.m_model(mod_index - 1).m_chain(j).chain_type() == pinang::water ||
                    pdb1.m_model(mod_index - 1).m_chain(j).chain_type() == pinang::ion)
                    continue;
                for (int n = 0; n < pdb1.m_model(mod_index - 1).m_chain(j).m_chain_length(); n++) {
                    r1 = pdb1.m_model(mod_index - 1).m_chain(j).m_residue(n);
                    cti1 = r1.chain_type();
                    std::vector<pinang::Residue> resi_group1;
                    if (cti1 == pinang::protein) {
                        resi_group1.push_back(r1);
                    }
                    else if (cti1 == pinang::DNA || cti1 == pinang::RNA || cti1 == pinang::na) {
                        pinang::Atom atmp = r1.m_P();
                        pinang::Residue rtmp_P, rtmp_S, rtmp_B;
                        rtmp_P.set_resid_index(r1.resid_index());
                        rtmp_S.set_resid_index(r1.resid_index());
                        rtmp_B.set_resid_index(r1.resid_index());
                        rtmp_P.set_resid_name("P");
                        rtmp_S.set_resid_name("S");
                        rtmp_B.set_resid_name(r1.resid_name());
                        for (int v = 0; v < r1.m_residue_size(); v++) {
                            atmp = r1.m_atom(v);
                            std::string aname = atmp.atom_name();
                            if (aname == "P  " || aname == "OP1"
                                || aname == "OP2" || aname == "O5'")
                                rtmp_P.add_atom(atmp);
                            else if (aname[2] == '\'') {
                                rtmp_S.add_atom(atmp);
                            } else {
                                rtmp_B.add_atom(atmp);
                            }
                        }
                        if (n != 0)
                            resi_group1.push_back(rtmp_P);
                        resi_group1.push_back(rtmp_S);
                        resi_group1.push_back(rtmp_B);
                    }

                    // -------------------- Calculating distances --------------
                    for (int q = 0; q < int(resi_group1.size()); q++) {
                        pinang::Residue rr1 = resi_group1[q];
                        double dist_min = pinang::resid_min_distance(r0, rr1);

                        if (cti1 == pinang::protein && (min_dist_W_PRO == 0 || min_dist_W_PRO > dist_min))
                            min_dist_W_PRO = dist_min;

                        if (rr1.resid_name() == "P" && (min_dist_W_P == 0 || min_dist_W_P > dist_min))
                            min_dist_W_P = dist_min;
                        if (rr1.resid_name() == "S" && (min_dist_W_S == 0 || min_dist_W_S > dist_min))
                            min_dist_W_S = dist_min;

                        if (rr1.resid_name() == "A" && (min_dist_W_A == 0 || min_dist_W_A > dist_min))
                            min_dist_W_A = dist_min;
                        if (rr1.resid_name() == "T" && (min_dist_W_T == 0 || min_dist_W_T > dist_min))
                            min_dist_W_T = dist_min;
                        if (rr1.resid_name() == "G" && (min_dist_W_G == 0 || min_dist_W_G > dist_min))
                            min_dist_W_G = dist_min;
                        if (rr1.resid_name() == "C" && (min_dist_W_C == 0 || min_dist_W_C > dist_min))
                            min_dist_W_C = dist_min;
                    }
                }
            }
            if (min_dist_W_PRO < 100)
                out_file << " WAT_PAIR " << "PRO   "
                         << std::setw(6) << min_dist_W_PRO
                         << std::endl;
            if (min_dist_W_P < 100)
                out_file << " WAT_PAIR " << "P     "
                         << std::setw(6) << min_dist_W_P
                         << std::endl;
            if (min_dist_W_S < 100)
                out_file << " WAT_PAIR " << "S     "
                         << std::setw(6) << min_dist_W_S
                         << std::endl;
            if (min_dist_W_A < 100)
                out_file << " WAT_PAIR " << "A     "
                         << std::setw(6) << min_dist_W_A
                         << std::endl;
            if (min_dist_W_T < 100)
                out_file << " WAT_PAIR " << "T     "
                         << std::setw(6) << min_dist_W_T
                         << std::endl;
            if (min_dist_W_G < 100)
                out_file << " WAT_PAIR " << "G     "
                         << std::setw(6) << min_dist_W_G
                         << std::endl;
            if (min_dist_W_C < 100)
                out_file << " WAT_PAIR " << "C     "
                         << std::setw(6) << min_dist_W_C
                         << std::endl;
        }
    }


    return 0;
}

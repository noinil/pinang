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
    for (int i = 0; i < pdb1.m_model(mod_index - 1).m_model_size(); i++) {
        if (pdb1.m_model(mod_index - 1).m_chain(i).chain_type() != pinang::water)
            continue;
        for (int m = 0; m < pdb1.m_model(mod_index - 1).m_chain(i).m_chain_length(); m++) {
            r0 = pdb1.m_model(mod_index - 1).m_chain(i).m_residue(m);

            std::vector<double> min_dist_w_pro; // min distance water -- protein
            std::vector<pinang::Residue> spec;

            for (int j = 0; j < pdb1.m_model(mod_index - 1).m_model_size(); j++) {
                if (pdb1.m_model(mod_index - 1).m_chain(j).chain_type() != pinang::protein)
                    continue;
                double min_dist_W_PRO = 0;
                pinang::Residue special0;
                for (int n = 0; n < pdb1.m_model(mod_index - 1).m_chain(j).m_chain_length(); n++) {
                    r1 = pdb1.m_model(mod_index - 1).m_chain(j).m_residue(n);
                    double dist_min = pinang::resid_min_distance(r0, r1);
                    if (min_dist_W_PRO == 0 || min_dist_W_PRO > dist_min){
                        min_dist_W_PRO = dist_min;
                        special0 = r1;
                    }
                }
                if (min_dist_W_PRO < 10 && min_dist_W_PRO > 2){
                    min_dist_w_pro.push_back(min_dist_W_PRO);
                    spec.push_back(special0);
                    out_file << " WAT_PAIR " << r1.resid_name() << " "
                             << std::setw(6) << min_dist_W_PRO
                             << std::endl;
                }
            }

            int len1 = min_dist_w_pro.size();
            int len2 = spec.size();
            if (len1 != len2)
            {
                std::cout << " ==================== ERROR! "
                          << " Vector size mismatch! -------------------- "
                          << std::endl;
            }
            for (int s = 0; s < len1; s++) {
                double dist1 = min_dist_w_pro[s];
                if (dist1 > 5)
                    continue;
                pinang::Residue res1 = spec[s];
                for (int t = s + 1; t < len1; t++) {
                    double dist2 = min_dist_w_pro[t];
                    if (dist2 > 5) continue;
                    pinang::Residue res2 = spec[t];
                    double dist_PP = pinang::resid_min_distance(res1, res2);
                    out_file << " WAT_MED_PRO  " << res1.resid_name() << " "
                             << res2.resid_name() << " "
                             << std::setw(6) << dist_PP << std::endl;
                }
            }
            min_dist_w_pro.clear();
            spec.clear();
        }
    }

    return 0;
}

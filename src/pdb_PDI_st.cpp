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

    std::string out_name = "PDI_dist.dat";

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
                      << " -f some.pdb [-o PDI_dist.dat] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o PDI_dist.dat] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " ERROR: need parameter for option -f: " << std::endl
                  << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o PDI_dist.dat] [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    std::ofstream out_file(out_name.c_str());

    if (mod_flag != 1) {
        if (pdb1.get_n_models() == 1)
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
    for (int i = 0; i < pdb1.get_model(mod_index - 1).get_model_size(); i++) {
        c0 = c0 + pdb1.get_model(mod_index - 1).get_chain(i);
    }
    pinang::Residue r0, r1;
    pinang::ChainType cti0, cti1;
    for (int i = 0; i < pdb1.get_model(mod_index - 1).get_model_size(); i++) {
        for (int m = 0; m < pdb1.get_model(mod_index - 1).get_chain(i).get_chain_length(); m++) {
            r0 = pdb1.get_model(mod_index - 1).get_chain(i).get_residue(m);
            cti0 = r0.get_chain_type();
            if (cti0 == pinang::water || cti0 == pinang::ion)
                continue;
            std::vector<pinang::Atom> atom_group0;
            std::vector<pinang::Residue> resi_group0;
            if (cti0 == pinang::protein) {
                pinang::Atom atmp = r0.get_C_alpha();
                atom_group0.push_back(atmp);
                resi_group0.push_back(r0);
            }
            else if (cti0 == pinang::DNA || cti0 == pinang::RNA || cti0 == pinang::na) {
                pinang::Atom atmp = r0.get_P();
                atmp.set_resid_name("P");
                if (m != 0)
                    atom_group0.push_back(atmp);
                atmp = r0.get_S();
                atmp.set_resid_name("S");
                atom_group0.push_back(atmp);
                atmp = r0.get_B();
                atom_group0.push_back(atmp);

                pinang::Residue rtmp_P, rtmp_S, rtmp_B;
                rtmp_P.set_resid_index(r0.get_resid_index());
                rtmp_S.set_resid_index(r0.get_resid_index());
                rtmp_B.set_resid_index(r0.get_resid_index());
                for (int v = 0; v < r0.get_residue_size(); v++) {
                    atmp = r0.get_atom(v);
                    std::string aname = atmp.get_atom_name();
                    if (aname == "P  " || aname == "OP1"
                        || aname == "OP2")
                        rtmp_P.add_atom(atmp);
                    else if (aname[2] == '\'') {
                        rtmp_S.add_atom(atmp);
                    } else {
                        rtmp_B.add_atom(atmp);
                    }
                }
                if (m != 0)
                    resi_group0.push_back(rtmp_P);
                resi_group0.push_back(rtmp_S);
                resi_group0.push_back(rtmp_B);
            }
            for (int j = i + 1; j < pdb1.get_model(mod_index - 1).get_model_size(); j++) {
                for (int n = 0; n < pdb1.get_model(mod_index - 1).get_chain(j).get_chain_length(); n++) {
                    r1 = pdb1.get_model(mod_index - 1).get_chain(j).get_residue(n);
                    cti1 = r1.get_chain_type();
                    if (cti1 == pinang::water || cti1 == pinang::ion)
                        continue;
                    if ((cti0 == pinang::DNA || cti0 == pinang::na) &&
                        (cti1 == pinang::DNA || cti1 == pinang::na))
                        break;
                    std::vector<pinang::Atom> atom_group1;
                    std::vector<pinang::Residue> resi_group1;
                    if (cti1 == pinang::protein) {
                        pinang::Atom atmp = r1.get_C_alpha();
                        atom_group1.push_back(atmp);
                        resi_group1.push_back(r1);
                    }
                    else if (cti1 == pinang::DNA || cti1 == pinang::RNA || cti1 == pinang::na) {
                        pinang::Atom atmp = r1.get_P();
                        atmp.set_resid_name("P");
                        if (n != 0)
                            atom_group1.push_back(atmp);
                        atmp = r1.get_S();
                        atmp.set_resid_name("S");
                        atom_group1.push_back(atmp);
                        atmp = r1.get_B();
                        atmp.set_resid_name(r1.get_resid_name());
                        atom_group1.push_back(atmp);

                        pinang::Residue rtmp_P, rtmp_S, rtmp_B;
                        rtmp_P.set_resid_index(r1.get_resid_index());
                        rtmp_S.set_resid_index(r1.get_resid_index());
                        rtmp_B.set_resid_index(r1.get_resid_index());
                        for (int v = 0; v < r1.get_residue_size(); v++) {
                            atmp = r1.get_atom(v);
                            std::string aname = atmp.get_atom_name();
                            if (aname == "P  " || aname == "OP1"
                                || aname == "OP2")
                                rtmp_P.add_atom(atmp);
                            else if (aname[2] == '\'') {
                                rtmp_S.add_atom(atmp);
                            } else {
                                rtmp_B.add_atom(atmp);
                            }
                        }
                        // std::cout << rtmp_P << std::endl;
                        // std::cout << rtmp_S.get_residue_size() << std::endl;
                        // std::cout << rtmp_B << std::endl;
                        if (n != 0)
                            resi_group1.push_back(rtmp_P);
                        resi_group1.push_back(rtmp_S);
                        resi_group1.push_back(rtmp_B);
                    }

                    // -------------------- Calculating distances --------------
                    if (resi_group0.size() != atom_group0.size()) {
                        std::cout << "ERROR in getting atom group 0 and resi group 0"
                                  << std::endl;
                        return(1);
                    }
                    if (resi_group1.size() != atom_group1.size()) {
                        std::cout << "ERROR in getting atom group 1 and resi group 1"
                                  << std::endl;
                        return(1);
                    }
                    for (int p = 0; p < int(resi_group0.size()); p++) {
                        pinang::Atom a0 = atom_group0[p];
                        pinang::Residue rr0 = resi_group0[p];
                        for (int q = 0; q < int(resi_group1.size()); q++) {
                            pinang::Atom a1 = atom_group1[q];
                            pinang::Residue rr1 = resi_group1[q];
                            double dist_min = pinang::resid_min_distance(rr0, rr1);
                            double cut_off = 6.5;
                            if (a1.get_resid_name() == "P") cut_off = 5.5;
                            if (a1.get_resid_name() == "S") cut_off = 6.5;
                            if (a1.get_resid_name() == "A") cut_off = 6.6;
                            if (a1.get_resid_name() == "T") cut_off = 6.8;
                            if (a1.get_resid_name() == "G") cut_off = 6.7;
                            if (a1.get_resid_name() == "C") cut_off = 6.7;
                            if (dist_min < cut_off && dist_min > 0)
                            {
                                out_file << " RESID_PAIR " << std::setw(3) << i << " "
                                         << std::setw(6) << m << " "
                                         << std::setw(6) << r0.get_resid_name()
                                         << "   -- "
                                         << std::setw(3) << j << " "
                                         << std::setw(6) << n << " "
                                         << std::setw(6) << a1.get_resid_name()
                                         << "   dist_min = "
                                         << std::setw(6) << dist_min
                                         << std::endl;
                                double dist_Ca = pinang::atom_distance(a0, a1);
                                out_file << " CG_PAIR "
                                         << std::setw(5)<< a0.get_resid_name()
                                         << std::setw(5)<< a1.get_resid_name()
                                         << "  "
                                         << std::setw(6) << dist_Ca
                                         << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }


    return 0;
}

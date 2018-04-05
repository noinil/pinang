#include "PDB.hpp"

#include <iostream>
#include <iomanip>
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
            << "\n";
  std::cout << " ~          PINANG protein-DNA distance statistics        ~ "
            << "\n";
  std::cout << " ========================================================== "
            << "\n";

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
                  << "\n";
        exit(EXIT_SUCCESS);
        break;
      default: /* '?' */
        std::cout << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o PDI_dist.dat] [-m module] [-h]"
                  << "\n";
        exit(EXIT_FAILURE);
    }
  }

  if (!in_flag)
  {
    std::cout << " ERROR: need parameter for option -f: " << "\n"
              << " Usage: "
              << argv[0]
              << " -f some.pdb [-o PDI_dist.dat] [-m module] [-h]"
              << "\n";
    exit(EXIT_SUCCESS);
  }
  pinang::PDB pdb1(infilename);

  std::ofstream out_file(out_name.c_str());

  if (mod_flag != 1) {
    if (pdb1.get_size() == 1)
    {
      mod_index = 1;
    } else {
      // std::cout << " Please choose a MODULE: " ;
      // std::cin >> mod_index;
      std::cout << " Trying to use MODULE: 1 " ;
      mod_index = 1;
    }
  }

  std::cout << " Analyzing DNA curvature of MODULE " << mod_index
            << " of " << infilename  << " ... "
            << "\n";

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // main
  //  _ _                 _      _ _
  // ( | )_ __ ___   __ _(_)_ __( | )
  //  V V| '_ ` _ \ / _` | | '_ \V V
  //     | | | | | | (_| | | | | |
  //     |_| |_| |_|\__,_|_|_| |_|
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */
  pinang::Model mdl0 = pdb1.get_model(mod_index - 1);
  pinang::Chain c0;
  for (int i = 0; i < mdl0.get_size(); i++) {
    c0 = c0 + mdl0.get_chain(i);
  }
  pinang::Residue r0, r1;
  pinang::ChainType cti0, cti1;
  for (int i = 0; i < mdl0.get_size(); i++) {
    for (int m = 0; m < mdl0.get_chain(i).get_size(); m++) {
      r0 = mdl0.get_chain(i).get_residue(m);
      cti0 = r0.get_chain_type();
      if (cti0 == pinang::water || cti0 == pinang::ion)
        continue;
      std::vector<pinang::Atom> atom_group0;
      std::vector<pinang::Residue> resi_group0;
      if (cti0 == pinang::protein) {
        pinang::Atom atmp = r0.get_cg_C_alpha();
        atmp.set_residue_name(r0.get_residue_name());
        atom_group0.push_back(atmp);
        resi_group0.push_back(r0);
      }
      else if (cti0 == pinang::DNA || cti0 == pinang::RNA || cti0 == pinang::na) {
        pinang::Atom atmp = r0.get_cg_P();
        atmp.set_residue_name("P");
        if (m != 0)
          atom_group0.push_back(atmp);
        atmp = r0.get_cg_S();
        atmp.set_residue_name("S");
        atom_group0.push_back(atmp);
        atmp = r0.get_cg_B();
        atom_group0.push_back(atmp);

        pinang::Residue rtmp_P, rtmp_S, rtmp_B;
        rtmp_P.set_residue_serial(r0.get_residue_serial());
        rtmp_S.set_residue_serial(r0.get_residue_serial());
        rtmp_B.set_residue_serial(r0.get_residue_serial());
        for (int v = 0; v < r0.get_size(); v++) {
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
      for (int j = i + 1; j < mdl0.get_size(); j++) {
        for (int n = 0; n < mdl0.get_chain(j).get_size(); n++) {
          r1 = mdl0.get_chain(j).get_residue(n);
          cti1 = r1.get_chain_type();
          if (cti1 == pinang::water || cti1 == pinang::ion)
            continue;
          if ((cti0 == pinang::DNA || cti0 == pinang::na) &&
              (cti1 == pinang::DNA || cti1 == pinang::na))
            break;
          if (cti0 == pinang::protein && cti1 == pinang::protein )
            break;
          std::vector<pinang::Atom> atom_group1;
          std::vector<pinang::Residue> resi_group1;
          if (cti1 == pinang::protein) {
            pinang::Atom atmp = r1.get_cg_C_alpha();
            atom_group1.push_back(atmp);
            resi_group1.push_back(r1);
          }
          else if (cti1 == pinang::DNA || cti1 == pinang::RNA || cti1 == pinang::na) {
            pinang::Atom atmp = r1.get_cg_P();
            atmp.set_residue_name("P");
            if (n != 0)
              atom_group1.push_back(atmp);
            atmp = r1.get_cg_S();
            atmp.set_residue_name("S");
            atom_group1.push_back(atmp);
            atmp = r1.get_cg_B();
            atmp.set_residue_name(r1.get_residue_name());
            atom_group1.push_back(atmp);

            pinang::Residue rtmp_P, rtmp_S, rtmp_B;
            rtmp_P.set_residue_serial(r1.get_residue_serial());
            rtmp_S.set_residue_serial(r1.get_residue_serial());
            rtmp_B.set_residue_serial(r1.get_residue_serial());
            for (int v = 0; v < r1.get_size(); v++) {
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
            // std::cout << rtmp_P << "\n";
            // std::cout << rtmp_S.get_size() << "\n";
            // std::cout << rtmp_B << "\n";
            if (n != 0)
              resi_group1.push_back(rtmp_P);
            resi_group1.push_back(rtmp_S);
            resi_group1.push_back(rtmp_B);
          }

          // -------------------- Calculating distances --------------
          if (resi_group0.size() != atom_group0.size()) {
            std::cout << "ERROR in getting atom group 0 and resi group 0"
                      << "\n";
            return(1);
          }
          if (resi_group1.size() != atom_group1.size()) {
            std::cout << "ERROR in getting atom group 1 and resi group 1"
                      << "\n";
            return(1);
          }
          for (int p = 0; p < int(resi_group0.size()); p++) {
            pinang::Atom a0 = atom_group0[p];
            pinang::Residue rr0 = resi_group0[p];
            for (int q = 0; q < int(resi_group1.size()); q++) {
              pinang::Atom a1 = atom_group1[q];
              pinang::Residue rr1 = resi_group1[q];
              double dist_min = pinang::residue_min_distance(rr0, rr1);
              double cut_off = 8.0;
              // if (a1.get_residue_name() == "P") cut_off = 5.5;
              // if (a1.get_residue_name() == "S") cut_off = 6.5;
              // if (a1.get_residue_name() == "A") cut_off = 6.6;
              // if (a1.get_residue_name() == "T") cut_off = 6.8;
              // if (a1.get_residue_name() == "G") cut_off = 6.7;
              // if (a1.get_residue_name() == "C") cut_off = 6.7;
              if (dist_min < cut_off && dist_min > 0)
              {
                out_file << " RESID_PAIR " << std::setw(3) << i << " "
                         << std::setw(6) << m << " "
                         << std::setw(6) << r0.get_residue_name()
                         << "   -- "
                         << std::setw(3) << j << " "
                         << std::setw(6) << n << " "
                         << std::setw(6) << r1.get_residue_name()
                         << "   dist_min = "
                         << std::setw(6) << dist_min
                         << "\n";
                double dist_Ca = pinang::atom_distance(a0, a1);
                out_file << " CG_PAIR "
                         << std::setw(5)<< a0.get_residue_name()
                         << std::setw(5)<< a1.get_residue_name()
                         << "  "
                         << std::setw(6) << dist_Ca
                         << "\n";
              }
            }
          }
        }
      }
    }
  }


  return 0;
}

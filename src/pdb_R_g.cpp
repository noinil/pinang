#include "PDB.h"
#include "constants.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <boost/algorithm/string.hpp>

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
  int inp_flag = 0;

  std::string infilename = "some.pdb";
  std::string out_name = "rg.dat";
  std::string inp_name = "protein.inp";

  while ((opt = getopt(argc, argv, "o:x:b:g:i:m:f:h")) != -1) {
    switch (opt) {
      case 'o':
        out_name = optarg;
        break;
      case 'i':
        inp_name = optarg;
        inp_flag = 1;
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
                  << " -f some.pdb [-i protein_seq.inp] [-o rg.dat] \n"
                  << " [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
        break;
      default: /* '?' */
        std::cout << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-i protein_seq.inp] [-o rg.dat] \n"
                  << " [-m module] [-h]"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
  }

  if (!in_flag)
  {
    std::cout << " ERROR: need parameter for option -f and -i: "
              << std::endl << " Usage: " << argv[0]
              << " -f some.pdb [-i protein_seq.inp] [-o rg.dat] \n"
              << " [-m module] [-h]"
              << std::endl;
    exit(EXIT_SUCCESS);
  }

  pinang::PDB pdb1(infilename);

  std::ifstream inp_file(inp_name.c_str());
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

  std::cout << " Analyzing Radius of Gyration (R_g) of MODULE " << mod_index
            << " of " << infilename  << " ... " << std::endl
            << std::endl;

  /* ============================================================
  //      _                   _
  //     (_)_ __  _ __  _   _| |_
  //     | | '_ \| '_ \| | | | __|
  //     | | | | | |_) | |_| | |_
  //     |_|_| |_| .__/ \__,_|\__|
  //              |_|
  // ============================================================
  */
  std::cout << " 1. Read in structures ..." << std::endl;
  int i = 0;
  int j = 0;

  std::vector<char> chain_id;
  std::string inp_line;
  std::string tmp_str;
  while (inp_file.good()) {
    std::getline(inp_file, inp_line);

    if (inp_file.fail())
    {
      break;
    }

    tmp_str = inp_line.substr(0,6);

    if (tmp_str == "CHAIN:")
    {
      std::string tmp_s;
      std::istringstream tmp_sstr;
      inp_line.erase(0,6);
      std::vector<std::string> strs;
      boost::split(strs, inp_line, boost::is_any_of(","));
      for (i = 0; i < int(strs.size()); i++) {
        tmp_s = strs[i];
        char tmp_c = 0;
        tmp_sstr.str(strs[i]);
        tmp_sstr >> tmp_c;
        chain_id.push_back(tmp_c);
        tmp_sstr.clear();
      }
    }
  }
  if (chain_id.size() == 0)
  {
    std::cout << " WARNING: STRAND info not found!" << std::endl;
    chain_id.push_back('A');
    // exit(EXIT_FAILURE);
  }

  inp_file.close();

  std::vector<pinang::Vec3d> c_alpha_coors;
  pinang::Model mdl0 = pdb1.m_model(mod_index-1);
  pinang::Chain chain_tmp;
  int mdl_size = mdl0.m_model_size();
  for (i = 0; i < int(chain_id.size()); ++i)
    for (j = 0; j < mdl_size; ++j)
      if (mdl0.m_chain(j).chain_ID() == chain_id[i]) {
        chain_tmp = mdl0.m_chain(j);
        int len1 = chain_tmp.m_chain_length();
        for (int k = 1; k < len1; k++) {
          c_alpha_coors.push_back(chain_tmp.m_residue(k).m_C_alpha().coordinates());
        }
        chain_tmp.reset();
      }

  std::cout << " ... Calculating ... " << std::endl;

  int len_ca_coors = c_alpha_coors.size();
  std::cout << "number of C_alpha: " << len_ca_coors << std::endl;

  pinang::Vec3d com_ca(0,0,0);
  for (i = 0; i < len_ca_coors; ++i) {
    com_ca = com_ca + c_alpha_coors[i];
  }
  com_ca = com_ca * (1.0 / len_ca_coors);

  double dist_square_sum = 0;
  for (i = 0; i < len_ca_coors; ++i) {
    double dist = pinang::vec_distance(c_alpha_coors[i], com_ca);
    std::cout << i << "  " << dist << std::endl;
    dist_square_sum += dist * dist;
  }

  double rg = sqrt(dist_square_sum / len_ca_coors);
  std::cout << rg << std::endl;

  std::cout << " ... done." << std::endl;


  // ----------------------------------------------------------------------
  out_file.close();

  return 0;
}


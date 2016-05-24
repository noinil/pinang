#include "PDB.hpp"
#include "constants.hpp"
#include "vec3d.hpp"

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace std;

void print_usage(char* s);

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
  int pdbin_flag = 0;
  int inp_flag = 0;
  int out_flag = 0;

  std::string infilename = "some.pdb";
  std::string out_name = "rg.dat";
  std::string inp_name = "group.inp";

  while ((opt = getopt(argc, argv, "o:i:m:f:h")) != -1) {
    switch (opt) {
      case 'o':
        out_name = optarg;
        out_flag = 1;
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
        pdbin_flag = 1;
        break;
      case 'h':
        print_usage(argv[0]);
        break;
      default: /* '?' */
        print_usage(argv[0]);
    }
  }

  if (!pdbin_flag)
  {
    std::cout << " ERROR: need parameter for option -f and -i: ";
    print_usage(argv[0]);
  }

  pinang::PDB pdb1(infilename);


  if (mod_flag != 1) {
    if (pdb1.get_n_models() == 1)
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
  std::vector<char> chain_id;

  int i = 0;
  int j = 0;


  if (inp_flag == 1) {
    std::ifstream inp_file(inp_name.c_str());
    std::string inp_line;
    std::string tmp_str;
    while (inp_file.good()) {
      std::getline(inp_file, inp_line);
      if (inp_file.fail())
        break;
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
    inp_file.close();
  }
  if (chain_id.size() == 0)
  {
    std::cout << " WARNING: STRAND info not found!" << std::endl;
    chain_id.push_back('A');
    // exit(EXIT_FAILURE);
  }


  std::vector<pinang::Vec3d> c_alpha_coors;
  pinang::Model mdl0 = pdb1.get_model(mod_index-1);
  pinang::Chain chain_tmp;
  int mdl_size = mdl0.get_model_size();
  for (i = 0; i < int(chain_id.size()); ++i)
    for (j = 0; j < mdl_size; ++j)
      if (mdl0.get_chain(j).get_chain_ID() == chain_id[i]) {
        chain_tmp = mdl0.get_chain(j);
        int len1 = chain_tmp.get_chain_length();
        for (int k = 1; k < len1; k++) {
          c_alpha_coors.push_back(chain_tmp.get_residue(k).get_cg_C_alpha().get_coordinate());
        }
        chain_tmp.reset();
      }

  std::cout << " ... Calculating ... " << std::endl;

  int len_ca_coors = c_alpha_coors.size();
  std::cout << "number of C_alpha: " << len_ca_coors << std::endl;

  pinang::Vec3d com_ca(0,0,0);
  pinang::Vec3d com_all(0,0,0);
  for (i = 0; i < len_ca_coors; ++i) {
    com_all = com_all + c_alpha_coors[i];
  }
  com_ca = com_all * (1.0 / len_ca_coors);

  double dist_square_sum = 0;
  for (i = 0; i < len_ca_coors; ++i) {
    pinang::Vec3d v_dist = c_alpha_coors[i] - com_ca;
    double sqr_dist = v_dist.squared_norm();
    dist_square_sum += sqr_dist;
  }

  double rg = sqrt(dist_square_sum / len_ca_coors);
  std::cout << rg << std::endl;

  std::cout << " ... done." << std::endl;


  // ----------------------------------------------------------------------
  if (out_flag == 1) {
    std::ofstream out_file(out_name.c_str());
    out_file << rg << std::endl;
    out_file.close();
  }


  return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: "
            << s
            << "\n\t -f some.pdb\n\t"
            << " [-i protein_seq.inp]"
            << "\n\t [-o rg.dat]\n\t"
            << " [-m module]\n\t [-h]"
            << std::endl;
  exit(EXIT_SUCCESS);
}

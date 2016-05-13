#include "PDB.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            << std::endl;
  std::cout << " ~            PINANG PDB CG positions output              ~ "
            << std::endl;
  std::cout << " ========================================================== "
            << std::endl;

  int opt, mod_index = 0;
  int mod_flag = 0;
  int in_flag = 0;

  std::string infilename = "some.pdb";
  std::string pos_name = "out.pos";
  std::string top_name = "out.top";

  while ((opt = getopt(argc, argv, "p:t:m:f:h")) != -1) {
    switch (opt) {
      case 'p':
        pos_name = optarg;
        break;
      case 't':
        top_name = optarg;
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
        print_usage(argv[0]);
        break;
      default: /* '?' */
        print_usage(argv[0]);
    }
  }

  if (!in_flag)
  {
    std::cout << " ERROR: need parameter for option -f: " << std::endl;
    print_usage(argv[0]);
  }
  pinang::PDB pdb1(infilename);

  std::ofstream pos_file(pos_name.c_str());
  std::ofstream top_file(top_name.c_str());

  if (mod_flag != 1) {
    if (pdb1.get_n_models() == 1)
    {
      mod_index = 1;
    } else {
      std::cout << " Please choose a MODULE: " ;
      std::cin >> mod_index;
    }
  }

  std::cout << " Extracting C-alpha coordinates of MODULE " << mod_index
            << " of " << infilename
            << " to " << pos_name << " ... "
            << std::endl;
  pos_file << "# CG positions for PDB " << infilename << std::endl;
  pdb1.get_model(mod_index - 1).output_cg_pos(pos_file);
  std::cout << " Done! " << std::endl
            << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            << std::endl;

  std::cout << " Extracting topological information of MODULE "
            << mod_index
            << " of " << infilename
            << " to " << top_name << " ... "
            << std::endl;
  top_file << "# CG topology for PDB " << infilename << std::endl;
  pdb1.get_model(mod_index - 1).output_top_mass(top_file);
  pdb1.get_model(mod_index - 1).output_top_bond(top_file);
  pdb1.get_model(mod_index - 1).output_top_angle(top_file);
  pdb1.get_model(mod_index - 1).output_top_dihedral(top_file);
  pdb1.get_model(mod_index - 1).output_top_nonbonded(top_file);
  // std::cout << "everything OK by now......" << std::endl;
  std::cout << " Done! " << std::endl
            << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            << std::endl;

  return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: "
            << s
            << "\n\t -f some.pdb\n\t [-p out.pos]\n\t"
            << " [-t out.top]\n\t [-m module]\n\t [-h]"
            << std::endl;
  exit(EXIT_SUCCESS);
}

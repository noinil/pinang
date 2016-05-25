/*!
@file pdb_cg_top.cpp
@brief Generate .top and .pos files from PDB structures.

Read PDB file, extract information for molecules, and generate .top / .pos files
for additional analysis.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-24 18:11
@copyright GNU Public License V3.0
*/


#include <fstream>
#include <unistd.h>
#include "PDB.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
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
    std::cout << " ERROR: need parameter for option -f: " << "\n";
    print_usage(argv[0]);
  }
  pinang::PDB pdb1(infilename);

  std::ofstream pos_file(pos_name.c_str());
  std::ofstream top_file(top_name.c_str());

  if (mod_flag != 1) {
    if (pdb1.get_size() == 1)
    {
      mod_index = 0;
    } else {
      std::cout << " Please choose a MODULE: " ;
      std::cin >> mod_index;
      --mod_index;
    }
  }

  pos_file << "# CG positions for PDB " << infilename << "\n";
  pdb1.get_model(mod_index).output_cg_pos(pos_file);

  top_file << "# CG topology for PDB " << infilename << "\n";
  pdb1.get_model(mod_index).output_top_mass(top_file);
  pdb1.get_model(mod_index).output_top_bond(top_file);
  pdb1.get_model(mod_index).output_top_angle(top_file);
  pdb1.get_model(mod_index).output_top_dihedral(top_file);
  pdb1.get_model(mod_index).output_top_nonbonded(top_file);

  pos_file.close();
  top_file.close();

  return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: "
            << s
            << "\n\t -f some.pdb\n\t [-p out.pos]\n\t"
            << " [-t out.top]\n\t [-m module]\n\t [-h]"
            << "\n";
  exit(EXIT_SUCCESS);
}

/*!
@file pdb_cg_top.cpp
@brief Generate .psf and .pdb files from PDB structures.

Read PDB file, extract information for molecules, and generate .psf / .pdb files
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
  std::string crd_name = "cg.pdb";
  std::string top_name = "cg.psf";
  std::string parm_name = "cg.ffp";

  while ((opt = getopt(argc, argv, "p:t:c:m:f:h")) != -1) {
    switch (opt) {
      case 'p':
        parm_name = optarg;
        break;
      case 'c':
        crd_name = optarg;
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

  std::ofstream crd_file(crd_name.c_str());
  std::ofstream top_file(top_name.c_str());
  std::ofstream parm_file(parm_name.c_str());

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

  pinang::Model m0 = pdb1.get_model(mod_index);
  m0.output_cg_crd(crd_file);
  crd_file << "END" << std::endl;

  m0.output_top_mass(top_file);
  m0.output_top_bond(top_file);
  m0.output_top_angle(top_file);
  m0.output_top_dihedral(top_file);

  m0.output_ffparm_bond(parm_file);
  m0.output_ffparm_angle(parm_file);
  m0.output_ffparm_dihedral(parm_file);
  m0.output_ffparm_nonbonded(parm_file);

  crd_file.close();
  top_file.close();
  parm_file.close();

  return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: "
            << s
            << "\n\t -f some.pdb\n\t [-p out.crd]\n\t"
            << " [-t out.top]\n\t [-m module]\n\t [-h]"
            << "\n";
  exit(EXIT_SUCCESS);
}

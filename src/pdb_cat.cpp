/*!
************************************************************
@file pdb_cat.cpp
@brief Re-output PDB structures.

Extract model from PDB, and reoutput the ATOM records to new PDB file.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-24 18:10
@copyright GNU Public License V3.0
************************************************************
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
  std::string outfilename = "out.pdb";

  while ((opt = getopt(argc, argv, "o:m:f:h")) != -1) {
    switch (opt) {
      case 'o':
        outfilename = optarg;
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
    print_usage(argv[0]);
  }
  pinang::PDB pdb1(infilename);

  std::ofstream ofile(outfilename.c_str());
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

  ofile << pdb1.get_model(mod_index);
  ofile << "END" << "\n";
  ofile.close();

  return 0;
}

void print_usage(char* s)
{
  std::cout << "  Usage: "
            << s
            << "\n\t -f some.pdb\n\t [-o output.pdb]\n\t [-m module]\n\t [-h]"
            << "\n";
  std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            << "\n";
  exit(EXIT_SUCCESS);
}

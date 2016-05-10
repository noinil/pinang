#include "PDB.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
  std::cout << " ~           PINANG sequence print            ~ " << std::endl;
  std::cout << " ============================================== " << std::endl;

  std::string infilename = "some.pdb";
  std::string out_name = "seq.fasta";

  int opt;
  int in_flag = 0;
  int out_flag = 0;

  while ((opt = getopt(argc, argv, "o:f:h")) != -1) {
    switch (opt) {
      case 'o':
        out_name = optarg;
        out_flag = 1;
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


  std::cout << " Sequence of PDB "
            << pdb1.get_pdb_name()
            << " :"
            << std::endl;
  std::cout << " Total number of chains: "
            << pdb1.get_model(0).get_model_size()
            << std::endl;
  std::cout << std::endl;
  std::cout << " 1-char-aa-name : ----------------------------- " << std::endl;
  pdb1.print_sequence(1);
  std::cout << " ---------------------------------------------- " << std::endl;
  std::cout << " 3-char-aa-name : ----------------------------- " << std::endl;
  pdb1.print_sequence(3);
  std::cout << " ---------------------------------------------- " << std::endl;
  if (out_flag) {
    std::ofstream out_file(out_name.c_str());
    pdb1.output_fasta(out_file);
  }

  return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: "
            << s
            << "\n\t -f some.pdb\n\t [-o seq.fasta]\n\t [-h]"
            << std::endl;
  exit(EXIT_SUCCESS);
}

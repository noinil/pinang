/*!
  @file pdb_get_sequence.cpp
  @brief Output sequence of protein/NA from PDB.

  Read PDB, extract sequences of protein/DNA/RNA, and output sequence to screen or
  file.  Both 3-char name and 1-char name are output.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 18:06
  @copyright GNU Public License V3.0
*/

#include <fstream>
#include <unistd.h>
#include "PDB.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  string infilename = "some.pdb";
  string out_name = "seq.fasta";

  int opt;
  int in_flag = 0;
  int out_flag = 0;

  while ((opt = getopt(argc, argv, "of:h")) != -1) {
    switch (opt) {
      case 'o':
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
    cout << " ERROR: need parameter for option -f: " << "\n";
    print_usage(argv[0]);
  }
  pinang::PDB pdb1(infilename);


  cout << " Sequence of PDB " << pdb1.get_pdb_name() << " :" << "\n";
  cout << " Total number of chains: " << pdb1.get_model(0).get_size() << "\n";
  cout << "\n";
  cout << " 1-char-aa-name : ----------------------------- " << "\n";
  pdb1.output_sequence(1);
  cout << " ---------------------------------------------- " << "\n";
  cout << " 3-char-aa-name : ----------------------------- " << "\n";
  pdb1.output_sequence(3);
  cout << " ---------------------------------------------- " << "\n";
  cout << "\n";

  if (out_flag) {
    out_name = infilename.substr(0, infilename.size()-4);
    out_name += ".fasta";
    ofstream out_file(out_name.c_str());
    pdb1.output_sequence_fasta(out_file);
    out_file.close();
  }

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
       << s
       << "\n\t -f some.pdb\n\t [-o seq.fasta]\n\t [-h]"
       << "\n";
  exit(EXIT_SUCCESS);
}

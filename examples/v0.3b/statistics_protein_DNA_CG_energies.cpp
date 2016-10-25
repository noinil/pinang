/*!
  @file statistics_protein_DNA_CG_energies.cpp
  @brief Statistics of protein-DNA sequence-specific interaction energies.

  Calculate protein-DNA sequence-specific interaction energies.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-06-16 14:33
  @copyright GNU Public License V3.0
*/

#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include "geometry.hpp"
#include "topology.hpp"
#include "ff_protein_DNA_specific.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  int opt, mod_index = 0;
  int in_flag = 0;

  string basefilename = "";
  string infilename = "some.psf";
  string outfilename;

  while ((opt = getopt(argc, argv, "o:f:h")) != -1) {
    switch (opt) {
      case 'o':
        outfilename = optarg;
        break;
      case 'f':
        infilename = optarg;
        in_flag = 1;
        basefilename = infilename.substr(0, infilename.size()-4);
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

  // -------------------- READ IN --------------------
  pinang::Topology top(infilename);

  string crd_name = basefilename + ".crd";
  pinang::Conformation conf(crd_name);


  string ffp_name = basefilename + ".ffp";
  outfilename = basefilename + "_pro_DNA_energy.dat";
  ofstream outfile(outfilename.c_str());


  // 0 energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " ============================================================" << "\n";
  pinang::FFProteinDNASpecific ff_ss(ffp_name);
  cout << " ============================================================" << "\n";
  double total_energy_0 = 0;
  total_energy_0 += ff_ss.compute_energy_protein_DNA_specific(top, conf);
  cout << "Total energy (native): " << total_energy_0 << "\n";

}

void print_usage(char* s)
{
  cout << " Usage: "
       << s
       << "\n\t -f some.psf\n\t"
       << " [-h]"
       << "\n";
  exit(EXIT_SUCCESS);
}



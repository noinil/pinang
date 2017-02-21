/*!
  @file cafedcd_Ep.cpp
  @brief Recalculate potential energies from CafeMol DCD file.

  Read DCD (CafeMol) file, recalculate potential energies based on force field parameters.

  @author Cheng Tan (noinil@gmail.com)
  @date 2017-02-21 11:39
  @copyright GNU Public License V3.0
*/

#include "read_cafemol_dcd.hpp"
#include "topology.hpp"
#include "geometry.hpp"
#include "ff_protein_DNA_specific.hpp"

#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <boost/algorithm/string.hpp>

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  int opt;
  int out_flag = 0;

  string dcd_name = "please_provide_name.dcd";
  string top_name = "please_provide_name.psf";
  string ffp_name = "please_provide_name.ffp";
  string ene_name = "please_provide_name.dat";
  string basefilename = "";

  while ((opt = getopt(argc, argv, "f:s:p:o:h")) != -1) {
    switch (opt) {
      case 'f':
        dcd_name = optarg;
        break;
      case 's':
        top_name = optarg;
        basefilename = top_name.substr(0, top_name.size()-4);
        break;
      case 'p':
        ffp_name = optarg;
        break;
      case 'o':
        ene_name = optarg;
        out_flag = 1;
        break;
      case 'h':
        print_usage(argv[0]);
        break;
      default: /* '?' */
        print_usage(argv[0]);
    }
  }

  // ------------------------------ prepare files ------------------------------
  if (out_flag == 0) {
    ene_name = basefilename + "_Ep.dat";
  }
  ifstream dcd_file(dcd_name.c_str(), ifstream::binary);
  ofstream ene_file(ene_name.c_str());
  pinang::Topology top(top_name);
  pinang::Conformation conf_tmp;

  // ------------------------------ Reading DCD --------------------------------
  vector<pinang::Conformation> conformations;
  pinang::read_cafemol_dcd(dcd_file, conformations);
  int nframe = conformations.size();
  if (nframe == 0)
  {
    cout << " ERROR: Empty DCD file!  Please check! " << "\n";
    return 1;
  }
  if (top.get_size() != conformations[0].get_size())
  {
    cout << " ERROR: Particle number don't match in top and dcd! "
         << " Please check! " << "\n";
    return 1;
  }

  // ------------------------------ Calculating energies --------------------------
  double total_energy_0 = 0;
  double ene_pdss = 0;
  double ene_ele = 0;

  pinang::FFProteinDNASpecific ff_ss(ffp_name);

  cout << " Calculating energies from dcd file : " << dcd_name << " ... " << endl;
  for (int i= 0; i < nframe; ++i) {
    conf_tmp = conformations[i];
    // ------------------------------ PDSS ------------------------------
    ene_pdss = ff_ss.compute_energy_protein_DNA_specific(top, conf_tmp);
    cout << "Protein-DNA sequence-specific energy: " << ene_pdss << "\n";
  }

  dcd_file.close();
  ene_file.close();

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
       << s
       << " -f xxx.dcd -s xxx.psf -p xxx.ffp [-o xxx_Ep.dat] [-h]"
       << endl;
  exit(EXIT_SUCCESS);
}

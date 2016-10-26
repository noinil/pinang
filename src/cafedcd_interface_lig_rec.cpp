/*!
  @file cafedcd_interface_lig_rec.cpp
  @brief Find interface (contact) between rec and lig, from MD trajectory (dcd file).

  Read DCD (CafeMol) file, extract contacts between receptor and ligand.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-10-26 20:50
  @copyright GNU Public License V3.0
*/

#include "read_cafemol_dcd.hpp"
#include "topology.hpp"
#include "group.hpp"
#include "selection.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <boost/algorithm/string.hpp>

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  double contact_cutoff = 10.0;
  double com_cutoff = 0.0;

  int opt;

  string dcd_name = "please_provide_name.dcd";
  string top_name = "please_provide_name.psf";
  string inp_name = "please_provide_name.in";
  string dat_name = "please_provide_name.dat";

  while ((opt = getopt(argc, argv, "f:c:C:s:i:o:h")) != -1) {
    switch (opt) {
      case 'f':
        dcd_name = optarg;
        break;
      case 's':
        top_name = optarg;
        break;
      case 'c':
        contact_cutoff = atof(optarg);
        break;
      case 'C':
        com_cutoff = atof(optarg);
        break;
      case 'i':
        inp_name = optarg;
        break;
      case 'o':
        dat_name = optarg;
        break;
      case 'h':
        print_usage(argv[0]);
        break;
      default: /* '?' */
        print_usage(argv[0]);
    }
  }

  // ------------------------------ prepare files ------------------------------
  ifstream dcd_file(dcd_name.c_str(), ifstream::binary);
  ofstream dat_file(dat_name.c_str());
  pinang::Topology top(top_name);

  // ------------------------------ get selections -----------------------------
  pinang::Selection sel_lig(inp_name, "LIG");
  pinang::Selection sel_rec(inp_name, "REC");

  cout << " Number of particles in GROUP LIG: " << sel_lig.get_size() << "\n";
  cout << " Number of particles in GROUP REC: " << sel_rec.get_size() << "\n";

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

  // ------------------------------ Estimate size of rec/lig -------------------
  pinang::Group grp_lig_0(conformations[0], sel_lig);
  pinang::Group grp_rec_0(conformations[0], sel_rec);
  double rg_lig_0 = pinang::get_radius_of_gyration(grp_lig_0);
  double rg_rec_0 = pinang::get_radius_of_gyration(grp_rec_0);
  com_cutoff = 2 * (rg_lig_0 + rg_rec_0);

  // ------------------------------ Calculating interface ----------------------
  pinang::Vec3d centroid_lig;
  pinang::Vec3d centroid_rec;
  pinang::Vec3d coor_lig;
  pinang::Vec3d coor_rec;
  double d_tmp;
  vector<int> lig_resid_flag;
  vector<int> rec_resid_flag;
  for (int i= 0; i < nframe; ++i) {
    pinang::Group grp_lig(conformations[i], sel_lig);
    pinang::Group grp_rec(conformations[i], sel_rec);
    centroid_lig = grp_lig.get_centroid();
    centroid_rec = grp_rec.get_centroid();
    d_tmp = pinang::vec_distance(centroid_rec, centroid_lig);
    if (d_tmp > com_cutoff) {
      continue;
    }

    for (int j = 0; j < sel_rec.get_size(); ++j) {
      rec_resid_flag.push_back(0);
      lig_resid_flag.push_back(0);
    }

    for (int j = 0; j < sel_rec.get_size(); ++j) {
      coor_rec = conformations[i].get_coordinate(sel_rec.get_selection(j));
      for (int k = 0; k < sel_lig.get_size(); ++k) {
        coor_lig = conformations[i].get_coordinate(sel_lig.get_selection(k));
        d_tmp = pinang::vec_distance(coor_rec, coor_lig);
        if (d_tmp < contact_cutoff) {
          rec_resid_flag[j] = 1;
          lig_resid_flag[k] = 1;
        }
      }
    }
    dat_file << "STEP > " << setw(8) << i << " \n";
    dat_file << " | LIG > ";
    for (int j = 0; j < sel_rec.get_size(); ++j) {
      if (lig_resid_flag[j] > 0)
        dat_file << " " << setw(5) << sel_lig.get_selection(j) + 1;
    }
    dat_file << "\n | REC > ";
    for (int k = 0; k < sel_rec.get_size(); ++k) {
      if (rec_resid_flag[k] > 0)
        dat_file << " " << setw(5) << sel_rec.get_selection(k) + 1;
    }
    dat_file << "\nTER >" << "\n";

    rec_resid_flag.clear();
    lig_resid_flag.clear();
  }

  dcd_file.close();
  dat_file.close();

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
       << s
       << " -f xxx.dcd -s xxx.psf -i xxx.in \n"
       << "[-C com_cutoff(real)] [-c contact_cutoff(real)] [-o distance.dat] [-h]"
       << "\n";
  cout << " Input file example: \n"
       << " ~~~~~~~~~~~~~~~~~~~~ \n REC: 1 to 100 \n LIG: 2 to 50, 55 to 106 \n"
       << endl;
  exit(EXIT_SUCCESS);
}

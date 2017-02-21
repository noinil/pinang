/*!
  @file cafedcd_distance_lig_rec.cpp
  @brief Get minimum distance from COM of lig to rec from MD trajectory (dcd file).

  Read DCD (CafeMol) file, extract minimum distance from COM of lig to rec.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-10-25 17:00
  @copyright GNU Public License V3.0
*/

#include "read_cafemol_dcd.hpp"
#include "topology.hpp"
#include "group.hpp"

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

  string dcd_name = "please_provide_name.dcd";
  string top_name = "please_provide_name.psf";
  string inp_name = "please_provide_name.in";
  string dis_name = "please_provide_name.dat";

  while ((opt = getopt(argc, argv, "f:s:i:o:h")) != -1) {
    switch (opt) {
      case 'f':
        dcd_name = optarg;
        break;
      case 's':
        top_name = optarg;
        break;
      case 'i':
        inp_name = optarg;
        break;
      case 'o':
        dis_name = optarg;
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
  ofstream dis_file(dis_name.c_str());
  pinang::Topology top(top_name);

  // ------------------------------ get selections -----------------------------
  pinang::Selection sel_lig(inp_name, "LIG");
  pinang::Selection sel_rec(inp_name, "REC");

  cout << " Number of particles in GROUP LIG: " << sel_lig.get_size() << "\n";
  cout << " Number of particles in GROUP REC: " << sel_rec.get_size() << "\n";

  vector<double> masses_lig;
  for (int i = 0; i < sel_lig.get_size(); ++i) {
    masses_lig.push_back(top.get_particle(sel_lig.get_selection(i)).get_mass());
  }

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

  // ------------------------------ Calculating distance ----------------------
  double dist;
  double d_tmp;
  pinang::Vec3d com_lig;
  pinang::Vec3d coor_rec;
  cout << " Calculating distance_min from LIG(COM) to REC : ..." << endl;
  for (int i= 0; i < nframe; ++i) {
    pinang::Group grp_lig(conformations[i], sel_lig);
    com_lig = pinang::get_center_of_mass(grp_lig, masses_lig);

    dist = -1.0;
    d_tmp = 0.0;
    for (int j = 0; j < sel_rec.get_size(); ++j) {
      coor_rec = conformations[i].get_coordinate(sel_rec.get_selection(j));
      d_tmp = pinang::vec_distance(com_lig, coor_rec);
      if (dist < 0 || dist > d_tmp) dist = d_tmp;
    }

    dis_file << setw(6) << i
             << "   " << setw(8) << dist
             << "\n"; // Output the distance!
  }
  cout << " Done! " << "\n";

  dcd_file.close();
  dis_file.close();

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
            << s
            << " -f xxx.dcd -s xxx.psf -i xxx.in [-o distance.dat] [-h]"
            << "\n";
  cout << " Input file example: \n"
       << " ~~~~~~~~~~~~~~~~~~~~ \n REC: 1 to 100 \n LIG: 2 to 50, 55 to 106 \n"
       << endl;
  exit(EXIT_SUCCESS);
}

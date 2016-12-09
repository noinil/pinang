/*!
  @file cafedcd_angle.cpp
  @brief Get angle between two vectors, calculating from MD trajectory (dcd file).

  Read DCD (CafeMol) file, extract angles formed between two vectors.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-12-08 15:52
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
  pinang::Selection sel_vecA_1(inp_name, "VA1");
  pinang::Selection sel_vecA_2(inp_name, "VA2");
  pinang::Selection sel_vecB_1(inp_name, "VB1");
  pinang::Selection sel_vecB_2(inp_name, "VB2");

  cout << " Number of particles in GROUP VEC_A_1: " << sel_vecA_1.get_size() << "\n";
  cout << " Number of particles in GROUP VEC_A_2: " << sel_vecA_2.get_size() << "\n";
  cout << " Number of particles in GROUP VEC_B_1: " << sel_vecB_1.get_size() << "\n";
  cout << " Number of particles in GROUP VEC_B_2: " << sel_vecB_2.get_size() << "\n";

  vector<double> masses_vecA_1;
  vector<double> masses_vecA_2;
  vector<double> masses_vecB_1;
  vector<double> masses_vecB_2;
  for (int i = 0; i < sel_vecA_1.get_size(); ++i) 
    masses_vecA_1.push_back(top.get_particle(sel_vecA_1.get_selection(i)).get_mass());
  for (int i = 0; i < sel_vecA_2.get_size(); ++i) 
    masses_vecA_2.push_back(top.get_particle(sel_vecA_2.get_selection(i)).get_mass());
  for (int i = 0; i < sel_vecB_1.get_size(); ++i) 
    masses_vecB_1.push_back(top.get_particle(sel_vecB_1.get_selection(i)).get_mass());
  for (int i = 0; i < sel_vecB_2.get_size(); ++i) 
    masses_vecB_2.push_back(top.get_particle(sel_vecB_2.get_selection(i)).get_mass());

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

  // ------------------------------ Calculating angle ----------------------
  double angle;
  pinang::Vec3d com_vecA_1;
  pinang::Vec3d com_vecA_2;
  pinang::Vec3d com_vecB_1;
  pinang::Vec3d com_vecB_2;
  pinang::Vec3d vecA;
  pinang::Vec3d vecB;
  cout << " Calculating angle : ..." << endl;
  for (int i= 0; i < nframe; ++i) {
    pinang::Group grp_vecA_1(conformations[i], sel_vecA_1);
    pinang::Group grp_vecA_2(conformations[i], sel_vecA_2);
    pinang::Group grp_vecB_1(conformations[i], sel_vecB_1);
    pinang::Group grp_vecB_2(conformations[i], sel_vecB_2);
    com_vecA_1 = pinang::get_center_of_mass(grp_vecA_1, masses_vecA_1);
    com_vecA_2 = pinang::get_center_of_mass(grp_vecA_2, masses_vecA_2);
    com_vecB_1 = pinang::get_center_of_mass(grp_vecB_1, masses_vecB_1);
    com_vecB_2 = pinang::get_center_of_mass(grp_vecB_2, masses_vecB_2);

    vecA = com_vecA_2 - com_vecA_1;
    vecB = com_vecB_2 - com_vecB_1;

    angle = vec_angle_deg(vecA, vecB);

    dis_file << setw(6) << i
             << "   " << setw(8) << angle
             << "\n"; // Output the angle!
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
            << " -f xxx.dcd -s xxx.psf -i xxx.in [-o angle.dat] [-h]"
            << "\n";
  cout << " Input file example (vec1 = VA2 - VA1; vec2 = VB2 - VB1): \n"
       << " ~~~~~~~~~~~~~~~~~~~~ \n VA1: 1 to 2 \n VA2: 3 to 50, 55 to 66 \n"
       << " VB1: 101 to 108 \n VB2: 120 to 160 \n"
       << endl;
  exit(EXIT_SUCCESS);
}

/*!
  @file statistics_protein_DNA_CG_mutation.cpp
  @brief Statistics of protein-DNA interacting quantities with mutation.

  Calculate protein-DNA interacting pairwise distances, angles, etc...

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-06-16 14:33
  @copyright GNU Public License V3.0
*/

#include <fstream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include "geometry.hpp"
#include "topology.hpp"

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

  vector<vector<int>> interaction_pairs;
  vector<int> pair(2, int());

  string stat_name = basefilename + ".stat";
  ifstream stat_file(stat_name.c_str());
  outfilename = basefilename + "_pro_DNA_local.dat";
  ofstream outfile(outfilename.c_str());

  string stat_line;
  string tmp_s;
  int tmp_i1, tmp_i2, tmp_i3, tmp_i4, tmp_i5, tmp_i6;
  string::size_type m;
  while (stat_file.good()) {
    getline(stat_file, stat_line);
    if (stat_file.fail())
      break;

    istringstream tmp_sstr;
    tmp_sstr.str(stat_line);
    m = stat_line.find("DB");
    if (m == 69) {
      for (int i = 0; i < 9; ++i)
        tmp_sstr >> tmp_s;
      tmp_sstr >> tmp_i1;
      for (int i = 0; i < 7; ++i)
        tmp_sstr >> tmp_s;
      tmp_sstr >> tmp_i2;
      pair[0] = tmp_i1 - 1;
      pair[1] = tmp_i2 - 1;
      if (top.get_particle(tmp_i2-1).get_atom_name() != "DB  ") {
        cout << " Wrong interaction pair for DB!!! \n";
        exit(EXIT_SUCCESS);
      }
      interaction_pairs.push_back(pair);
    }
  }

  // -------------------- Calculating... --------------------
  pinang::Vec3d tmp_c_CA, tmp_c_B0, tmp_c_S0, tmp_c_B5, tmp_c_B3, tmp_c_new_DB;
  pinang::Vec3d vec_BS, vec_BC, vec_53;
  double tmp_dist_0;
  double tmp_angle_0, tmp_angle_5, tmp_angle_3, tmp_angle_53;
  for (int i = 0; i < interaction_pairs.size(); ++i) {
    tmp_i1 = interaction_pairs[i][0];
    tmp_i2 = interaction_pairs[i][1];
    tmp_i3 = tmp_i2 - 1;
    tmp_i4 = tmp_i2 - 3;
    tmp_i5 = tmp_i2 + 3;

    if ( (tmp_i4 < 0 ||
          top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID())
         && (tmp_i5 >= top.get_size() ||
             top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()))
      continue;

    tmp_c_CA = conf.get_coordinate(tmp_i1);  // Coor of CA
    tmp_c_B0 = conf.get_coordinate(tmp_i2);  // Coor of B
    tmp_c_S0 = conf.get_coordinate(tmp_i3);  // Coor of S connected to B

    vec_BS = tmp_c_S0 - tmp_c_B0;
    vec_BC = tmp_c_CA - tmp_c_B0;

    tmp_dist_0 = vec_BC.norm();
    tmp_angle_0 = pinang::vec_angle_deg(vec_BS, vec_BC);

    if (tmp_i4 < 0
        || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      tmp_c_B5 = conf.get_coordinate(tmp_i2);
      // tmp_angle_5 = -999;
    } else {
      tmp_c_B5 = conf.get_coordinate(tmp_i4);  // Coor of 5' B
      // tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_B0, vec_BC);
    }
    if (tmp_i5 >= top.get_size()
        || top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      tmp_c_B3 = conf.get_coordinate(tmp_i2);
      // tmp_angle_3 = -999;
    } else {
      tmp_c_B3 = conf.get_coordinate(tmp_i5);  // Coor of 3' B
      // tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_B0, vec_BC);
    }
    vec_53 = tmp_c_B3 - tmp_c_B5;
    tmp_angle_53 = pinang::vec_angle_deg(vec_BC, vec_53);

    outfile << setw(4) << tmp_i1 + 1 << "  "
            << setw(4) << tmp_i2 + 1 << "   "
            << setiosflags(ios_base::fixed) << setprecision(6)
            << setw(10) << tmp_dist_0 << "  "
            << setiosflags(ios_base::fixed) << setprecision(3)
            << setw(9) << tmp_angle_0 << "   "
            << setw(9) << tmp_angle_53 << "\n";
  }

  stat_file.close();
  outfile.close();

  return 0;
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

/*!
  @file statistics_protein_DNA_CG_quantities.cpp
  @brief Statistics of protein-DNA interacting quantities.

  Calculate protein-DNA interacting pairwise distances, angles, etc...

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-06-16 14:33
  @copyright GNU Public License V3.0
*/

#include <fstream>
#include <sstream>
#include <unistd.h>
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
  outfilename = basefilename + "_pro_DNA_stat.dat";
  ofstream outfile(outfilename.c_str());

  string stat_line;
  string tmp_s;
  int tmp_i1, tmp_i2, tmp_i3, tmp_i4, tmp_i5;
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
  pinang::Vec3d tmp_c_CA, tmp_c_B0, tmp_c_S0, tmp_c_B5, tmp_c_B3, tmp_c_new_CA;
  pinang::Vec3d tmp_c_CA_N, tmp_c_CA_C;
  pinang::Vec3d tmp_c_new_CA_N, tmp_c_new_CA_C;
  double tmp_dist_0;
  double tmp_angle_0, tmp_angle_53, tmp_angle_NC;
  for (int i = 0; i < interaction_pairs.size(); ++i) {
    tmp_i1 = interaction_pairs[i][0];
    tmp_i2 = interaction_pairs[i][1];
    outfile << "------------------------------------------------------------ \n";
    outfile << "CONTACT PAIR: " << tmp_i1 + 1 << " - " << tmp_i2 + 1 << "\n";

    // ---------- Base -- CA distance ----------
    tmp_c_CA = conf.get_coordinate(tmp_i1);  // Coor of CA
    tmp_c_B0 = conf.get_coordinate(tmp_i2);  // Coor of B
    tmp_dist_0 = pinang::vec_distance(tmp_c_CA, tmp_c_B0);
    outfile << " Distance DB0-CA : " << tmp_dist_0 << "\n";

    int k = tmp_i1 - 1;
    int calpha_term_N = 0, calpha_term_C = 0;
    if (k < 0 || top.get_particle(k).get_chain_ID() != top.get_particle(tmp_i1).get_chain_ID()) {
      tmp_c_CA_N = tmp_c_CA;
      calpha_term_N = 1;
    } else {
      calpha_term_N = 0;
      tmp_c_CA_N = conf.get_coordinate(k);  // Coor of N' CA
    }
    k = i + 1;
    if (k >= top.get_size() || top.get_particle(k).get_chain_ID() != top.get_particle(tmp_i1).get_chain_ID()) {
      tmp_c_CA_C = tmp_c_CA;
      calpha_term_C = 1;
    } else {
      calpha_term_C = 0;
      tmp_c_CA_C = conf.get_coordinate(k);  // Coor of C' CA
    }
    if (calpha_term_N * calpha_term_C > 0) {
      std::cout << " Single Residue Protein Chain!!! WTF!!! \n";
      exit(EXIT_SUCCESS);
    }
    tmp_angle_NC = pinang::vec_angle_deg(tmp_c_CA_C - tmp_c_CA_N, tmp_c_B0 - tmp_c_CA);
    outfile << " Angle NCA-CCA--CA-DB0 : " << tmp_angle_NC << "\n";

    // ---------- Sugar -- Base -- CA angle ----------
    tmp_i3 = tmp_i2 - 1;
    if (top.get_particle(tmp_i3).get_atom_name() != "DS  ") {
      cout << " Wrong interaction pair for DB!!! \n";
      exit(EXIT_SUCCESS);
    }
    tmp_c_S0 = conf.get_coordinate(tmp_i3);  // Coor of S connected to B 
    tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    outfile << " Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";

    int term_5 = 0, term_3 = 0;
    // ---------- 5' Base ----------
    tmp_i4 = tmp_i2 - 3;
    if (tmp_i4 < 0 || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      term_5 = 1;
    } else {
      if (top.get_particle(tmp_i4).get_atom_name() != "DB  ") {
        cout << " Wrong interaction pair for 5'DB!!! \n";
        exit(EXIT_SUCCESS);
      }
      tmp_c_B5 = conf.get_coordinate(tmp_i4);  // Coor of 5' B
    }
    // ---------- 3' Base ----------
    tmp_i4 = tmp_i2 + 3;
    if (tmp_i4 >= top.get_size() || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      term_3 = 1;
    } else {
      if (top.get_particle(tmp_i4).get_atom_name() != "DB  ") {
        cout << " Wrong interaction pair for 3'DB!!! \n";
        exit(EXIT_SUCCESS);
      }
      tmp_c_B3 = conf.get_coordinate(tmp_i4);  // Coor of 3' B
    }
    if (term_5 + term_3 == 0) {
      tmp_angle_53 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_B5, tmp_c_CA - tmp_c_B0);
      outfile << " Angle DB5'-DB3'--DB0-CA : " << tmp_angle_53 << "\n";
    } else {
      outfile << " Angle DB5'-DB3'--DB0-CA : NaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    }


    // ------------------------------ Shifting... ------------------------------
    tmp_i4 = tmp_i2 - 6;
    tmp_i5 = tmp_i2 + 3;
    if ( tmp_i4 < 0 || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID() || 
         tmp_i5 >= top.get_size() || top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " No 5' shifting... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else {
      vector<pinang::Vec3d> v5;
      for (int l = 0; l < 8; ++l) 
        v5.push_back(conf.get_coordinate(tmp_i4 - 1 + l));
      pinang::Group g5(v5);
      vector<pinang::Vec3d> v0;
      for (int l = 0; l < 8; ++l) 
        v0.push_back(conf.get_coordinate(tmp_i2 - 4 + l));
      pinang::Group g0(v0);
      // do the transform...
      pinang::Transform t;
      find_transform(g0, g5, t);
      tmp_c_new_CA = t.apply(tmp_c_CA);
      tmp_c_new_CA_C = t.apply(tmp_c_CA_C);
      tmp_c_new_CA_N = t.apply(tmp_c_CA_N);
      // Calculation....
      tmp_dist_0 = pinang::vec_distance(tmp_c_new_CA, tmp_c_B0);
      outfile << " 5'-Shifted Distance DB0-CA : " << tmp_dist_0 << "\n";

      tmp_angle_NC = pinang::vec_angle_deg(tmp_c_new_CA_C - tmp_c_new_CA_N, tmp_c_B0 - tmp_c_new_CA);
      outfile << " 5'-Shifted Angle NCA-CCA--CA-DB0 : " << tmp_angle_NC << "\n";

      tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_new_CA - tmp_c_B0);
      outfile << " 5'-Shifted Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";

      if (term_5 + term_3 == 0) {
        tmp_angle_53 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_B5, tmp_c_new_CA - tmp_c_B0);
        outfile << " 5'-Shifted Angle DB5'-DB3'--DB0-CA : " << tmp_angle_53 << "\n";
      } else {
        outfile << " 5'-Shifted DB5'-DB3'--DB0-CA : NaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
      }
    }
    tmp_i4 = tmp_i2 - 3;
    tmp_i5 = tmp_i2 + 6;
    if ( tmp_i4 < 0 || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID() || 
         tmp_i5 >= top.get_size() || top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " No 3' shifting... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else {
      vector<pinang::Vec3d> v3;
      for (int l = 0; l < 8; ++l) 
        v3.push_back(conf.get_coordinate(tmp_i2 - 1 + l));
      pinang::Group g3(v3);
      vector<pinang::Vec3d> v0;
      for (int l = 0; l < 8; ++l) 
        v0.push_back(conf.get_coordinate(tmp_i2 - 4 + l));
      pinang::Group g0(v0);
      // do the transform...
      pinang::Transform t;
      find_transform(g0, g3, t);
      tmp_c_new_CA = t.apply(tmp_c_CA);
      tmp_c_new_CA_C = t.apply(tmp_c_CA_C);
      tmp_c_new_CA_N = t.apply(tmp_c_CA_N);
      // Calculation....
      tmp_dist_0 = pinang::vec_distance(tmp_c_new_CA, tmp_c_B0);
      outfile << " 3'-Shifted Distance DB0-CA : " << tmp_dist_0 << "\n";

      tmp_angle_NC = pinang::vec_angle_deg(tmp_c_new_CA_C - tmp_c_new_CA_N, tmp_c_B0 - tmp_c_new_CA);
      outfile << " 3'-Shifted Angle NCA-CCA--CA-DB0 : " << tmp_angle_NC << "\n";

      tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_new_CA - tmp_c_B0);
      outfile << " 3'-Shifted Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";

      if (term_5 + term_3 == 0) {
        tmp_angle_53 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_B5, tmp_c_new_CA - tmp_c_B0);
        outfile << " 3'-Shifted Angle DB5'-DB3'--DB0-CA : " << tmp_angle_53 << "\n";
      } else {
        outfile << " 3'-Shifted DB5'-DB3'--DB0-CA : NaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
      }
    }

    outfile << "END PAIR \n";
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

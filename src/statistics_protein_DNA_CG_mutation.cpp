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
  outfilename = basefilename + "_pro_DNA_muta.dat";
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

  // ------ read in standard DNA cg conformation ------
  int i;
  pinang::Conformation std_cg_dna;
  std::vector<pinang::Vec3d> std_dna_coors(6, pinang::Vec3d());
  string std_dna_pdb_name = "/home/noinil/Learningspace/pdb_test/cgdna.pdb";
  pinang::Atom atom_tmp;
  int icount = 0;
  ifstream reffile(std_dna_pdb_name.c_str());
  while (reffile.good()) {
    reffile >> atom_tmp;
    if (reffile.fail())
      break;
    std_dna_coors[icount] = atom_tmp.get_coordinate();
    ++icount;
  }
  std_cg_dna.set_conformation(std_dna_coors);
  pinang::Selection sA, sC, sG, sT;
  std::vector<int> vsa = {0, 1, 2};
  std::vector<int> vst = {0, 1, 3};
  std::vector<int> vsg = {0, 1, 4};
  std::vector<int> vsc = {0, 1, 5};
  sA.set_selection(vsa);
  sT.set_selection(vst);
  sG.set_selection(vsg);
  sC.set_selection(vsc);
  pinang::Group std_cg_A(std_cg_dna, sA);
  pinang::Group std_cg_C(std_cg_dna, sC);
  pinang::Group std_cg_G(std_cg_dna, sG);
  pinang::Group std_cg_T(std_cg_dna, sT);

  // -------------------- Calculating... --------------------
  pinang::Vec3d tmp_c_CA, tmp_c_B0, tmp_c_S0, tmp_c_B5, tmp_c_B3, tmp_c_new_DB;
  double tmp_dist_0;
  double tmp_angle_0, tmp_angle_5, tmp_angle_3;
  for (int i = 0; i < interaction_pairs.size(); ++i) {
    tmp_i1 = interaction_pairs[i][0];
    tmp_i2 = interaction_pairs[i][1];
    outfile << "------------------------------------------------------------ \n";
    outfile << "CONTACT PAIR: " << tmp_i1 + 1
            << " - " << tmp_i2 + 1
            << "  |  " << top.get_particle(tmp_i1).get_residue_name()
            << " - " << top.get_particle(tmp_i2).get_residue_name()
            << "\n";

    // ---------- Base -- CA distance ----------
    tmp_c_CA = conf.get_coordinate(tmp_i1);  // Coor of CA
    tmp_c_B0 = conf.get_coordinate(tmp_i2);  // Coor of B
    tmp_dist_0 = pinang::vec_distance(tmp_c_CA, tmp_c_B0);
    outfile << " Distance DB0-CA : " << tmp_dist_0 << "\n";

    // ---------- Sugar -- Base -- CA angle ----------
    tmp_i3 = tmp_i2 - 1;
    if (top.get_particle(tmp_i3).get_atom_name() != "DS  ") {
      cout << " Wrong interaction pair for DB!!! \n";
      exit(EXIT_SUCCESS);
    }
    tmp_c_S0 = conf.get_coordinate(tmp_i3);  // Coor of S connected to B
    tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    outfile << " Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";

    // ---------- 5' Base ----------
    tmp_i4 = tmp_i2 - 3;
    if (tmp_i4 < 0 || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " Angle DB5'-DB0-CA : NaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else {
      if (top.get_particle(tmp_i4).get_atom_name() != "DB  ") {
        cout << " Wrong interaction pair for 5'DB!!! \n";
        exit(EXIT_SUCCESS);
      }
      tmp_c_B5 = conf.get_coordinate(tmp_i4);  // Coor of 5' B
      // ---------- 5' Base -- Base -- CA angle ----------
      tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
      outfile << " Angle DB5'-DB0-CA : " << tmp_angle_5 << "\n";
    }
    // ---------- 3' Base ----------
    tmp_i5 = tmp_i2 + 3;
    if (tmp_i5 >= top.get_size() || top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " Angle DB3'-DB0-CA : NaN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else {
      if (top.get_particle(tmp_i5).get_atom_name() != "DB  ") {
        cout << " Wrong interaction pair for 3'DB!!! \n";
        exit(EXIT_SUCCESS);
      }
      tmp_c_B3 = conf.get_coordinate(tmp_i5);  // Coor of 3' B
      // ---------- 3' Base -- Base -- CA angle ----------
      tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
      outfile << " Angle DB3'-DB0-CA : " << tmp_angle_3 << "\n";
    }

    // ------------------------------ Mutating... ------------------------------
    tmp_i6 = tmp_i2 - 2;
    pinang::Vec3d std_cg_Base;
    if (tmp_i6 < 0 || top.get_particle(tmp_i6).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " No mutation allowed! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else if (tmp_i4 < 0 || top.get_particle(tmp_i4).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " No mutation allowed! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else if (tmp_i5 >= top.get_size() || top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
      outfile << " No mutation allowed! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n";
    } else {
      vector<pinang::Vec3d> v0;
      v0.push_back(conf.get_coordinate(tmp_i6));
      v0.push_back(conf.get_coordinate(tmp_i3));
      v0.push_back(conf.get_coordinate(tmp_i2));
      pinang::Group g0(v0);
      // do the transform...
      pinang::Transform t;
      // ------ A mutant ------
      if (top.get_particle(tmp_i2).get_residue_name() != "DA ") {
        find_transform(std_cg_A, g0, t);
        std_cg_Base = std_cg_A.get_coordinate(2);
        tmp_c_new_DB = t.apply(std_cg_Base);
        tmp_dist_0 = pinang::vec_distance(tmp_c_CA, tmp_c_new_DB);
        outfile << " DA - Mutated Distance DB0-CA : " << tmp_dist_0 << "\n";
        tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DA - Mutated Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";
        tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DA - Mutated Angle DB5'-DB0-CA : " << tmp_angle_5 << "\n";
        tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DA - Mutated Angle DB3'-DB0-CA : " << tmp_angle_3 << "\n";
      }
      // ------ C mutant ------
      if (top.get_particle(tmp_i2).get_residue_name() != "DC ") {
        find_transform(std_cg_C, g0, t);
        std_cg_Base = std_cg_C.get_coordinate(2);
        tmp_c_new_DB = t.apply(std_cg_Base);
        tmp_dist_0 = pinang::vec_distance(tmp_c_CA, tmp_c_new_DB);
        outfile << " DC - Mutated Distance DB0-CA : " << tmp_dist_0 << "\n";
        tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DC - Mutated Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";
        tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DC - Mutated Angle DB5'-DB0-CA : " << tmp_angle_5 << "\n";
        tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DC - Mutated Angle DB3'-DB0-CA : " << tmp_angle_3 << "\n";
      }
      // ------ G mutant ------
      if (top.get_particle(tmp_i2).get_residue_name() != "DG ") {
        find_transform(std_cg_G, g0, t);
        std_cg_Base = std_cg_G.get_coordinate(2);
        tmp_c_new_DB = t.apply(std_cg_Base);
        tmp_dist_0 = pinang::vec_distance(tmp_c_CA, tmp_c_new_DB);
        outfile << " DG - Mutated Distance DB0-CA : " << tmp_dist_0 << "\n";
        tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DG - Mutated Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";
        tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DG - Mutated Angle DB5'-DB0-CA : " << tmp_angle_5 << "\n";
        tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DG - Mutated Angle DB3'-DB0-CA : " << tmp_angle_3 << "\n";
      }
      // ------ T mutant ------
      if (top.get_particle(tmp_i2).get_residue_name() != "DT ") {
        find_transform(std_cg_T, g0, t);
        std_cg_Base = std_cg_T.get_coordinate(2);
        tmp_c_new_DB = t.apply(std_cg_Base);
        tmp_dist_0 = pinang::vec_distance(tmp_c_CA, tmp_c_new_DB);
        outfile << " DT - Mutated Distance DB0-CA : " << tmp_dist_0 << "\n";
        tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DT - Mutated Angle DS0-DB0-CA : " << tmp_angle_0 << "\n";
        tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DT - Mutated Angle DB5'-DB0-CA : " << tmp_angle_5 << "\n";
        tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_new_DB, tmp_c_CA - tmp_c_new_DB);
        outfile << " DT - Mutated Angle DB3'-DB0-CA : " << tmp_angle_3 << "\n";
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

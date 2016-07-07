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

using namespace std;

const double d_theta0 = 20;
const double d_theta5 = 20;
const double d_theta3 = 20;

double factorize_theta(double theta, double theta_0, int i);
double pro_DNA_specific_energy(double r, double theta0, double theta5, double theta3,
                               double r0, double theta0_0, double theta5_0, double theta3_0);

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

  vector<int> dna_index;
  vector<vector<int>> interaction_pairs;
  vector<int> pair(2, int());
  vector<double> interaction_sigma;
  vector<double> interaction_angle_0;
  vector<double> interaction_angle_5;
  vector<double> interaction_angle_3;

  string ffp_name = basefilename + ".ffp";
  ifstream ffp_file(ffp_name.c_str());
  outfilename = basefilename + "_pro_DNA_energy.dat";
  ofstream outfile(outfilename.c_str());

  string ffp_line;
  string tmp_s;
  int tmp_i1, tmp_i2, tmp_i0;
  double tmp_dist, tmp_a0, tmp_a5, tmp_a3;
  string::size_type m;
  while (ffp_file.good()) {
    getline(ffp_file, ffp_line);
    if (ffp_file.fail())
      break;

    istringstream tmp_sstr;
    tmp_sstr.str(ffp_line);
    m = ffp_line.find("[ protein-DNA seq-specific ]");
    if (m != string::npos) {
      for (int i = 0; i < 4; ++i)
        tmp_sstr >> tmp_s;
      tmp_sstr >> tmp_i0;  // read in total number of pro-DNA interactions;
      getline(ffp_file, ffp_line);  // read in a comment line;
      for (int i = 0; i < tmp_i0; ++i) {  // read in all interaction details...;
        getline(ffp_file, ffp_line);
        istringstream tmp_sstr1;
        tmp_sstr1.str(ffp_line);
        tmp_sstr1 >> tmp_i1;
        tmp_sstr1 >> tmp_i2;
        pair[0] = tmp_i1 - 1;
        pair[1] = tmp_i2 - 1;
        if (top.get_particle(tmp_i2-1).get_atom_name() != "DB  ") {
          cout << " Wrong interaction pair for DB!!! \n";
          exit(EXIT_SUCCESS);
        }
        interaction_pairs.push_back(pair);
        dna_index.push_back(tmp_i2 - 1);
        tmp_sstr1 >> tmp_dist;
        tmp_sstr1 >> tmp_a0;
        tmp_sstr1 >> tmp_a5;
        tmp_sstr1 >> tmp_a3;
        interaction_sigma.push_back(tmp_dist);
        interaction_angle_0.push_back(tmp_a0);
        interaction_angle_5.push_back(tmp_a5);
        interaction_angle_3.push_back(tmp_a3);
      }
    }
  }


  // Choose the reference group...
  // ==============================
  vector<int> reference_g;
  vector<int> reference_tmp;
  int tmp_reference_length = 0;
  sort(dna_index.begin(), dna_index.end());
  auto last = std::unique(dna_index.begin(), dna_index.end());
  dna_index.erase(last, dna_index.end());
  char cid = '_';
  for (int i : dna_index) {
    if (top.get_particle(i).get_chain_ID() != cid) {
      if (reference_tmp.size() > tmp_reference_length) {
        reference_g = reference_tmp;
        tmp_reference_length = reference_g.size();
      }
      reference_tmp.clear();
      cid = top.get_particle(i).get_chain_ID();
    }
    reference_tmp.push_back(i);
  }
  if (reference_tmp.size() > tmp_reference_length) {
    reference_g = reference_tmp;
    tmp_reference_length = reference_g.size();
  }
  // ==============================

  // Preparing for shifting...
  vector<int> sel_ref1;
  vector<int> sel_ref2;
  int refg_s = reference_g.size();
  for (int i = reference_g[0] - 1; i <= reference_g[refg_s - 1] - 3; ++i) 
    sel_ref1.push_back(i);
  for (int i = reference_g[0] + 2; i <= reference_g[refg_s - 1]; ++i) 
    sel_ref2.push_back(i);
  pinang::Selection sref1(sel_ref1);
  pinang::Selection sref2(sel_ref2);
  pinang::Group gref1(conf, sref1);
  pinang::Group gref2(conf, sref2);
  pinang::Transform t;
  pinang::Vec3d tmp_c_CA0, tmp_c_CA, tmp_c_B0, tmp_c_S0, tmp_c_B5, tmp_c_B3;
  double tmp_distance, tmp_angle_0, tmp_angle_5, tmp_angle_3;

  // 0 energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " ============================================================" << "\n";
  double total_energy_0 = 0;
  double energy_tmp = 0;
  for (int i = 0; i < interaction_pairs.size(); ++i) {
    int proi = interaction_pairs[i][0];
    int dnai = interaction_pairs[i][1];
    double sigma = interaction_sigma[i];
    double angle_0 = interaction_angle_0[i];
    double angle_5 = interaction_angle_5[i];
    double angle_3 = interaction_angle_3[i];
    tmp_c_CA = conf.get_coordinate(proi);
    tmp_c_B0 = conf.get_coordinate(dnai);
    tmp_distance = pinang::vec_distance(tmp_c_CA, tmp_c_B0);
    tmp_c_S0 = conf.get_coordinate(dnai - 1);
    tmp_angle_0 = vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    // ---------- 5' Base ----------
    if (angle_5 < -600) {
      tmp_angle_5 = angle_5;
    } else {
      tmp_c_B5 = conf.get_coordinate(dnai - 3);
      tmp_angle_5 = vec_angle_deg(tmp_c_B5 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    }
    // ---------- 3' Base ----------
    if (angle_3 < -600) {
      tmp_angle_3 = angle_3;
    } else {
      tmp_c_B3 = conf.get_coordinate(dnai + 3);
      tmp_angle_3 = vec_angle_deg(tmp_c_B3 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    }
    energy_tmp = pro_DNA_specific_energy(tmp_distance, tmp_angle_0, tmp_angle_5, tmp_angle_3, sigma, angle_0, angle_5, angle_3);
    total_energy_0 += energy_tmp;
    cout << "   " << proi << " - " << dnai << "   " << energy_tmp << "\n";
  }
  cout << "Total energy (native): " << total_energy_0 << "\n";

  // 5' shifting... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " ============================================================" << "\n";
  double total_energy_5 = 0;
  pinang::find_transform(gref2, gref1, t);
  for (int i = 0; i < interaction_pairs.size(); ++i) {
    int proi = interaction_pairs[i][0];
    int dnai = interaction_pairs[i][1];
    double sigma = interaction_sigma[i];
    double angle_0 = interaction_angle_0[i];
    double angle_5 = interaction_angle_5[i];
    double angle_3 = interaction_angle_3[i];
    tmp_c_CA0 = conf.get_coordinate(proi);
    tmp_c_CA = t.apply(tmp_c_CA0);
    tmp_c_B0 = conf.get_coordinate(dnai);
    tmp_distance = pinang::vec_distance(tmp_c_CA, tmp_c_B0);
    tmp_c_S0 = conf.get_coordinate(dnai - 1);
    tmp_angle_0 = vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    // ---------- 5' Base ----------
    if (angle_5 < -600) {
      tmp_angle_5 = angle_5;
    } else {
      tmp_c_B5 = conf.get_coordinate(dnai - 3);
      tmp_angle_5 = vec_angle_deg(tmp_c_B5 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    }
    // ---------- 3' Base ----------
    if (angle_3 < -600) {
      tmp_angle_3 = angle_3;
    } else {
      tmp_c_B3 = conf.get_coordinate(dnai + 3);
      tmp_angle_3 = vec_angle_deg(tmp_c_B3 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    }
    energy_tmp = pro_DNA_specific_energy(tmp_distance, tmp_angle_0, tmp_angle_5, tmp_angle_3, sigma, angle_0, angle_5, angle_3);
    total_energy_5 += energy_tmp;
    cout << "   " << proi << " - " << dnai << "   " << energy_tmp << "\n";
  }
  cout << "Total energy (5' shifting): " << total_energy_5 << "\n";

  // 3' shifting... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " ============================================================" << "\n";
  double total_energy_3 = 0;
  pinang::find_transform(gref1, gref2, t);
  for (int i = 0; i < interaction_pairs.size(); ++i) {
    int proi = interaction_pairs[i][0];
    int dnai = interaction_pairs[i][1];
    double sigma = interaction_sigma[i];
    double angle_0 = interaction_angle_0[i];
    double angle_5 = interaction_angle_5[i];
    double angle_3 = interaction_angle_3[i];
    tmp_c_CA0 = conf.get_coordinate(proi);
    tmp_c_CA = t.apply(tmp_c_CA0);
    tmp_c_B0 = conf.get_coordinate(dnai);
    tmp_distance = pinang::vec_distance(tmp_c_CA, tmp_c_B0);
    tmp_c_S0 = conf.get_coordinate(dnai - 1);
    tmp_angle_0 = pinang::vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    // ---------- 5' Base ----------
    if (angle_5 < -600) {
      tmp_angle_5 = angle_5;
    } else {
      tmp_c_B5 = conf.get_coordinate(dnai - 3);
      tmp_angle_5 = pinang::vec_angle_deg(tmp_c_B5 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    }
    // ---------- 3' Base ----------
    if (angle_3 < -600) {
      tmp_angle_3 = angle_3;
    } else {
      tmp_c_B3 = conf.get_coordinate(dnai + 3);
      tmp_angle_3 = pinang::vec_angle_deg(tmp_c_B3 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
    }
    energy_tmp = pro_DNA_specific_energy(tmp_distance, tmp_angle_0, tmp_angle_5, tmp_angle_3, sigma, angle_0, angle_5, angle_3);
    total_energy_3 += energy_tmp;
    cout << "   " << proi << " - " << dnai << "   " << energy_tmp << "\n";
  }
  cout << "Total energy (3' shifting): " << total_energy_3 << "\n";

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


double pro_DNA_specific_energy(double r, double theta0, double theta5, double theta3,
                               double r0, double theta0_0, double theta5_0, double theta3_0)
{
  double ene = 0;
  double r0_over_r = r0 / r;
  double r0_over_r_2 = r0_over_r * r0_over_r;
  double r0_over_r_4 = r0_over_r_2 * r0_over_r_2;
  double r0_over_r_6 = r0_over_r_4 * r0_over_r_2;
  double r0_over_r_12 = r0_over_r_6 * r0_over_r_6;
  double r0_over_r_10 = r0_over_r_6 * r0_over_r_4;

  ene = factorize_theta(theta0, theta0_0, 0)
        * factorize_theta(theta5, theta5_0, 5)
        * factorize_theta(theta3, theta3_0, 3)
        * (5 * r0_over_r_12 - 6 * r0_over_r_10);

  if (ene > 0) {
    return 0;
  }
  return ene;
}

double factorize_theta(double theta, double theta_0, int i)
{
  double fct = 0;
  double theta_cutoff;
  double delta_theta = abs(theta - theta_0);
  if (i == 0) {
    theta_cutoff = d_theta0;
  } else if (i == 5) {
    theta_cutoff = d_theta5;
  } else if (i == 3) {
    theta_cutoff = d_theta3;
  } else {
    cout << " ERROR: Wrong factorization_theta code!..." << "\n";
    exit(EXIT_SUCCESS);
  }

  if (delta_theta < theta_cutoff) {
    fct = 1;
  } else if (delta_theta < theta_cutoff * 2) {
    double cos_theta = cos(delta_theta / 180.0 * 3.14159);
    fct = 1 - cos_theta * cos_theta;
  } else {
    fct = 0;
  }

  return fct;
}

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
  ffp_file.close();


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
  vector<int> sel_all;
  int refg_s = reference_g.size();
  for (int i = reference_g[0] - 1; i <= reference_g[refg_s - 1] - 3; ++i) 
    sel_ref1.push_back(i);
  for (int i = reference_g[0] + 2; i <= reference_g[refg_s - 1]; ++i) 
    sel_ref2.push_back(i);
  for (int i = 0; i < conf.get_size(); ++i) 
    sel_all.push_back(i);

  pinang::Selection sref1(sel_ref1);
  pinang::Selection sref2(sel_ref2);
  pinang::Selection sall(sel_all);
  pinang::Group gref1(conf, sref1);
  pinang::Group gref2(conf, sref2);
  pinang::Group group0(conf, sall);
  pinang::Transform t;
  pinang::Vec3d tmp_c_CA0, tmp_c_CA, tmp_c_B0, tmp_c_S0, tmp_c_B5, tmp_c_B3;
  double tmp_distance, tmp_angle_0, tmp_angle_5, tmp_angle_3;

  // 0 energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pinang::FFProteinDNASpecific ff_ss(ffp_name);
  cout << " ============================================================" << "\n";
  double total_energy_0 = 0;
  total_energy_0 += ff_ss.compute_energy_protein_DNA_specific(top, conf);
  cout << "Total energy (native): " << total_energy_0 << "\n";

  // 5' shifting... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // cout << " ============================================================" << "\n";
  // double total_energy_5 = 0;
  // pinang::find_transform(gref2, gref1, t);
  // pinang::Group conf_new5 = t.apply(group0);
  // total_energy_5 += ff_ss.compute_energy_protein_DNA_specific(top, conf_new5);
  // cout << "Total energy (5' shifting): " << total_energy_5 << "\n";

  // 3' shifting... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // cout << " ============================================================" << "\n";
  // double total_energy_3 = 0;
  // pinang::find_transform(gref1, gref2, t);
  // pinang::Group conf_new3 = t.apply(group0);
  // total_energy_3 += ff_ss.compute_energy_protein_DNA_specific(top, conf_new3);
  // cout << "Total energy (3' shifting): " << total_energy_3 << "\n";

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



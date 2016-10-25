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
#include "ff_protein_DNA_specific.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  int opt, mod_index = 0;
  int in_flag = 0;
  int seq_flag = 0;

  string basefilename = "";
  string infilename = "some.psf";
  string outfilename;
  string seqfilename;

  while ((opt = getopt(argc, argv, "o:f:h:p:")) != -1) {
    switch (opt) {
      case 'o':
        outfilename = optarg;
        break;
      case 'p':
        seqfilename = optarg;
        seq_flag = 1;
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
  if (!seq_flag)
  {
    cout << " ERROR: need seq file for option -p: " << "\n";
    print_usage(argv[0]);
  }

  // -------------------- READ IN --------------------
  pinang::Topology top0(infilename);

  string crd_name = basefilename + ".crd";
  pinang::Conformation conf0(crd_name);


  string ffp_name = basefilename + ".ffp";
  pinang::FFProteinDNASpecific ff_ss(ffp_name);

  vector<string> mutat_seq;
  vector<int> strand1 = {133, 136, 139, 142, 145};
  vector<int> strand2 = {116, 113, 110, 107, 104};

  ifstream seq_file(seqfilename.c_str());
  outfilename = basefilename + "_pro_DNA_muta_ene.dat";
  ofstream outfile(outfilename.c_str());

  string seq_line;
  string tmp_s;
  int tmp_i1, tmp_i2;
  int tmp_i3, tmp_i4;
  int tmp_i5, tmp_i6;
  string::size_type m;
  while (seq_file.good()) {
    getline(seq_file, seq_line);
    if (seq_file.fail())
      break;

    istringstream tmp_sstr;
    tmp_sstr.str(seq_line);
    if (seq_line.size() > 20) {
      continue;
    }
    tmp_sstr >> tmp_s;
    mutat_seq.push_back(tmp_s);
  }


  // ------ read in standard DNA cg conformation ------
  int i;
  pinang::Conformation std_cg_dna;
  std::vector<pinang::Vec3d> std_dna_coors(6, pinang::Vec3d());
  string std_dna_pdb_name = "/Users/noinil/Learningspace/pdb_test/cgdna.pdb";
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
  // 0 energy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  double total_energy_0 = 0;
  total_energy_0 = ff_ss.compute_energy_protein_DNA_specific(top0, conf0);
  outfile << " Native Total energy (native): " << total_energy_0 << "\n";

 for (int i = 0; i < mutat_seq.size(); ++i) {
    outfile << " SEQ: " << mutat_seq[i] << " : ";
    // cout << " SEQ: " << mutat_seq[i] << " : \n";
    pinang::Topology top = top0;
    pinang::Conformation conf = conf0;
    pinang::Vec3d tmp_c_new_B1, tmp_c_new_B2;
    for (int j = 0; j < 5; ++j) {
      string tmp_base_name;
      if (mutat_seq[i][j] == 'A') {
        tmp_base_name = "DA ";
      } else if (mutat_seq[i][j] == 'C') {
        tmp_base_name = "DC ";
      } else if (mutat_seq[i][j] == 'G') {
        tmp_base_name = "DG ";
      } else if (mutat_seq[i][j] == 'T') {
        tmp_base_name = "DT ";
      } 
      tmp_i1 = strand1[j] - 1;
      tmp_i2 = strand2[j] - 1;
      tmp_i3 = tmp_i1 - 1;
      tmp_i4 = tmp_i2 - 1;
      tmp_i5 = tmp_i1 - 2;
      tmp_i6 = tmp_i2 - 2;

      // ------------------------------ Mutating... ------------------------------
      pinang::Vec3d std_cg_Base1, std_cg_Base2;
      if (tmp_i5 < 0 || top.get_particle(tmp_i5).get_chain_ID() != top.get_particle(tmp_i1).get_chain_ID()) {
        outfile << " No mutation allowed! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ";
        break;;
      } else if (tmp_i6 < 0 || top.get_particle(tmp_i6).get_chain_ID() != top.get_particle(tmp_i2).get_chain_ID()) {
        outfile << " No mutation allowed! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ";
        break;;
      } else {
        vector<pinang::Vec3d> v1;
        v1.push_back(conf.get_coordinate(tmp_i5));
        v1.push_back(conf.get_coordinate(tmp_i3));
        v1.push_back(conf.get_coordinate(tmp_i1));
        pinang::Group g1(v1);
        vector<pinang::Vec3d> v2;
        v2.push_back(conf.get_coordinate(tmp_i6));
        v2.push_back(conf.get_coordinate(tmp_i4));
        v2.push_back(conf.get_coordinate(tmp_i2));
        pinang::Group g2(v2);
        pinang::Transform t1;
        pinang::Transform t2;
        if (top.get_particle(tmp_i1).get_residue_name() == tmp_base_name) {
          continue;
        }
        if (tmp_base_name == "DA ") {
          find_transform(std_cg_A, g1, t1);
          find_transform(std_cg_T, g2, t2);
          std_cg_Base1 = std_cg_A.get_coordinate(2);
          std_cg_Base2 = std_cg_T.get_coordinate(2);
          tmp_c_new_B1 = t1.apply(std_cg_Base1);
          conf.set_coordinate(tmp_i1, tmp_c_new_B1);
          top.get_particle(tmp_i1).set_residue_name("DA ");
          tmp_c_new_B2 = t2.apply(std_cg_Base2);
          conf.set_coordinate(tmp_i2, tmp_c_new_B2);
          top.get_particle(tmp_i2).set_residue_name("DT ");
        }
        if (tmp_base_name == "DC ") {
          find_transform(std_cg_C, g1, t1);
          find_transform(std_cg_G, g2, t2);
          std_cg_Base1 = std_cg_C.get_coordinate(2);
          std_cg_Base2 = std_cg_G.get_coordinate(2);
          tmp_c_new_B1 = t1.apply(std_cg_Base1);
          conf.set_coordinate(tmp_i1, tmp_c_new_B1);
          top.get_particle(tmp_i1).set_residue_name("DC ");
          tmp_c_new_B2 = t2.apply(std_cg_Base2);
          conf.set_coordinate(tmp_i2, tmp_c_new_B2);
          top.get_particle(tmp_i2).set_residue_name("DG ");
        }
        if (tmp_base_name == "DG ") {
          find_transform(std_cg_G, g1, t1);
          find_transform(std_cg_C, g2, t2);
          std_cg_Base1 = std_cg_G.get_coordinate(2);
          std_cg_Base2 = std_cg_C.get_coordinate(2);
          tmp_c_new_B1 = t1.apply(std_cg_Base1);
          conf.set_coordinate(tmp_i1, tmp_c_new_B1);
          top.get_particle(tmp_i1).set_residue_name("DG ");
          tmp_c_new_B2 = t2.apply(std_cg_Base2);
          conf.set_coordinate(tmp_i2, tmp_c_new_B2);
          top.get_particle(tmp_i2).set_residue_name("DC ");
        }
        if (tmp_base_name == "DT ") {
          find_transform(std_cg_T, g1, t1);
          find_transform(std_cg_A, g2, t2);
          std_cg_Base1 = std_cg_T.get_coordinate(2);
          std_cg_Base2 = std_cg_A.get_coordinate(2);
          tmp_c_new_B1 = t1.apply(std_cg_Base1);
          conf.set_coordinate(tmp_i1, tmp_c_new_B1);
          top.get_particle(tmp_i1).set_residue_name("DT ");
          tmp_c_new_B2 = t2.apply(std_cg_Base2);
          conf.set_coordinate(tmp_i2, tmp_c_new_B2);
          top.get_particle(tmp_i2).set_residue_name("DA ");
        }
      }
    }
    total_energy_0 = ff_ss.compute_energy_protein_DNA_specific(top, conf);
    outfile << total_energy_0 << " \n";
  }


  seq_file.close();
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

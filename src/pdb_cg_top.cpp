/*!
  @file pdb_cg_top.cpp
  @brief Generate .psf and .pdb files from PDB structures.

  Read PDB file, extract information for molecules, and generate .psf / .pdb files
  for additional analysis.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 18:11
  @copyright GNU Public License V3.0
*/


#include <fstream>
#include <unistd.h>
#include "PDB.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  int opt, mod_index = 0;
  int mod_flag = 0;
  int in_flag = 0;
  int crd_flag = 0;
  int pdb_flag = 0;
  int parm_flag = 0;
  int stat_flag = 0;

  string basefilename = "";
  string infilename = "some.pdb";

  while ((opt = getopt(argc, argv, "spcPm:f:h")) != -1) {
    switch (opt) {
      case 's':
        stat_flag = 1;
        break;
      case 'P':
        parm_flag = 1;
        break;
      case 'p':
        pdb_flag = 1;
        break;
      case 'c':
        crd_flag = 1;
        break;
      case 'm':
        mod_index = atoi(optarg);
        mod_flag = 1;
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
  pinang::PDB pdb1(infilename);


  if (mod_flag != 1) {
    if (pdb1.get_size() == 1)
    {
      mod_index = 0;
    } else {
      cout << " Please choose a MODULE: " ;
      cin >> mod_index;
      --mod_index;
    }
  }

  pinang::Model m0 = pdb1.get_model(mod_index);

  if (pdb_flag) {
    string pdb_name = basefilename + "_cg.pdb";
    ofstream pdb_file(pdb_name.c_str());
    m0.output_cg_pdb(pdb_file);
    pdb_file << "END" << endl;
    pdb_file.close();
  }

  if (crd_flag) {
    string crd_name = basefilename + "_cg.crd";
    ofstream crd_file(crd_name.c_str());
    m0.output_cg_crd(crd_file);
    crd_file << "END" << endl;
    crd_file.close();
  }

  if (parm_flag) {
    string parm_name = basefilename + "_cg.ffp";
    ofstream parm_file(parm_name.c_str());
    m0.output_ffparm_bond(parm_file);
    m0.output_ffparm_angle(parm_file);
    m0.output_ffparm_dihedral(parm_file);
    m0.output_ffparm_nonbonded(parm_file);
    parm_file.close();
  }

  string top_name = basefilename + "_cg.psf";
  ofstream top_file(top_name.c_str());
  top_file << "PSF CMAP " << endl << endl;
  top_file << "      3 !NTITLE" << endl;
  top_file << "REMARKS Created by pinang." << endl;
  top_file << "REMARKS Based on " << basefilename << endl;
  top_file << "REMARKS Enjoy " << endl << endl;
  m0.output_top_mass(top_file);
  m0.output_top_bond(top_file);
  m0.output_top_angle(top_file);
  m0.output_top_dihedral(top_file);
  top_file.close();

  if (stat_flag) {
    string stat_name = basefilename + "_cg.stat";
    ofstream stat_file(stat_name.c_str());
    m0.output_statistics_pro_DNA_contact_pairs(stat_file);
    stat_file << "END" << endl;
    stat_file.close();
  }

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
       << s
       << "\n\t -f some.pdb\n\t [-c (out_cg.crd)]\n\t"
       << " [-P (out_cg.ffp)]\n\t"
       << " [-p (out_cg.pdb)]\n\t [-s (out_cg.stat)]\n\t"
       << " [-m module]\n\t [-h]"
       << "\n";
  exit(EXIT_SUCCESS);
}

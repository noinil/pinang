/*!
  @file cafedcd_rmsd.cpp
  @brief Calculate RMSD from MD trajectory (dcd file).

  Read DCD (CafeMol) file, caalculate RMSD for a group of particles.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-10-25 20:48
  @copyright GNU Public License V3.0
*/

#include "read_cafemol_dcd.hpp"
#include "topology.hpp"
#include "geometry.hpp"

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
  int ref_flag = 0;

  string dcd_name = "please_provide_name.dcd";
  string top_name = "please_provide_name.psf";
  string inp_name = "please_provide_name.in";
  string ref_name = "please_provide_name.crd";
  string rmsd_name = "please_provide_name.dat";

  while ((opt = getopt(argc, argv, "f:s:i:o:r:h")) != -1) {
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
        rmsd_name = optarg;
        break;
      case 'r':
        ref_name = optarg;
        ref_flag = 1;
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
  ofstream rmsd_file(rmsd_name.c_str());
  pinang::Topology top(top_name);
  pinang::Conformation conf_ref;
  pinang::Conformation conf_obj;
  if (ref_flag == 1) {
    conf_ref = pinang::Conformation(ref_name);
  }

  // ------------------------------ get selections -----------------------------
  pinang::Selection sel_tran_ref(inp_name, "TRAN_REF");
  pinang::Selection sel_tran_obj(inp_name, "TRAN_OBJ");
  pinang::Selection sel_rmsd_ref(inp_name, "RMSD_REF");
  pinang::Selection sel_rmsd_obj(inp_name, "RMSD_OBJ");

  cout << " Number of particles in GROUP TRANSLATE_REFERENCE: " << sel_tran_ref.get_size() << "\n";
  cout << " Number of particles in GROUP TRANSLATE_OBJECT: " << sel_tran_obj.get_size() << "\n";
  cout << " Number of particles in GROUP RMSD_REFERENCE: " << sel_rmsd_ref.get_size() << "\n";
  cout << " Number of particles in GROUP RMSD_OBJECT: " << sel_rmsd_obj.get_size() << "\n";
  if (sel_tran_ref.get_size() != sel_tran_obj.get_size())
  {
    cout << " Error: Inconsistent number of particles in superimposition (RMSD). \n";
    print_usage(argv[0]);
  }
  if (sel_rmsd_ref.get_size() != sel_rmsd_obj.get_size())
  {
    cout << " Error: Inconsistent number of particles in RMSD calculation. \n";
    print_usage(argv[0]);
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

  // ------------------------------ Calculating rmsd --------------------------
  if (ref_flag == 0) {
    conf_ref = conformations[0];
  }
  pinang::Transform t;
  pinang::Group translated_obj;
  double rmsd;
  cout << " Calculating rmsd from dcd file : " << dcd_name << " ... " << endl;
  for (int i= 0; i < nframe; ++i) {
    conf_obj = conformations[i];
    pinang::Group grp_tran_ref(conf_ref, sel_tran_ref);
    pinang::Group grp_tran_obj(conf_obj, sel_tran_obj);
    pinang::Group grp_rmsd_ref(conf_ref, sel_rmsd_ref);
    pinang::Group grp_rmsd_obj(conf_obj, sel_rmsd_obj);

    pinang::find_transform(grp_tran_obj, grp_tran_ref, t);
    translated_obj = t.apply(grp_rmsd_obj);
    rmsd = pinang::get_rmsd(grp_rmsd_ref, translated_obj);
    rmsd_file << setw(6) << i
             << "   " << setw(8) << rmsd
             << "\n"; // Output the rmsdtance!
  }

  dcd_file.close();
  rmsd_file.close();

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
            << s
            << " -f xxx.dcd -s xxx.psf -i xxx.in [-r reference_struct.crd] [-o xxx_rmsd.dat] [-h]"
            << "\n";
  cout << " Input file example: \n"
       << " ~~~~~~~~~~~~~~~~~~~~ \n TRAN_REF: 1 to 100 \n TRAN_OBJ: 2 to 50, 55 to 106 \n"
       << " RMSD_REF: 1 to 200 \n RMSD_OBJ 2 to 201 \n ~~~~~~~~~~~~~~~~~~~~ "
       << endl;
  exit(EXIT_SUCCESS);
}

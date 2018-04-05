#include <fstream>
#include <sstream>
#include <unistd.h>
#include "geometry.hpp"
#include "topology.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  string reffilename = "dna0.pdb";
  string infilename = "all.stat";
  string outfilename = "translated.stat";
  ifstream infile(infilename.c_str());
  ifstream reffile(reffilename.c_str());
  ofstream outfile(outfilename.c_str());

  // -------------------- READ standard sugar coors -------
  pinang::Group std_sugar;
  std::vector<pinang::Vec3d> stdsugar_coors(6, pinang::Vec3d());
  pinang::Atom atom_tmp;
  int i, j, k;
  while (reffile.good()) {
    reffile >> atom_tmp;
    if (reffile.fail())
      break;
    if (atom_tmp.get_atom_name() == "C5' ") {
      stdsugar_coors[0] = atom_tmp.get_coordinate();
    } else if (atom_tmp.get_atom_name() == "C4' ") {
      stdsugar_coors[1] = atom_tmp.get_coordinate();
    } else if (atom_tmp.get_atom_name() == "O4' ") {
      stdsugar_coors[2] = atom_tmp.get_coordinate();
    } else if (atom_tmp.get_atom_name() == "C3' ") {
      stdsugar_coors[3] = atom_tmp.get_coordinate();
    } else if (atom_tmp.get_atom_name() == "C2' ") {
      stdsugar_coors[4] = atom_tmp.get_coordinate();
    } else if (atom_tmp.get_atom_name() == "C1' ") {
      stdsugar_coors[5] = atom_tmp.get_coordinate();
    }
  }
  std_sugar.set_conformation(stdsugar_coors);

  // -------------------- READ IN --------------------
  pinang::Group tmp_sugar;
  pinang::Atom atmp;
  pinang::Vec3d vcoor, vnew;
  std::vector<pinang::Vec3d> tmpsugar_coors(6, pinang::Vec3d());
  std::vector<pinang::Atom> tmp_atoms;
  string in_line;
  while (infile.good()) {
    getline(infile, in_line);
    istringstream tmp_sstr;
    tmp_sstr.str(in_line);
    tmp_sstr >> atom_tmp;
    if (infile.fail())
      break;

    if (atom_tmp.get_record_name() == "REMARK") {
      outfile << in_line << std::endl;
      for (int i = 0; i < 6; ++i) {
        infile >> atom_tmp;
        if (atom_tmp.get_atom_name() == "C5' ") {
          tmpsugar_coors[0] = atom_tmp.get_coordinate();
        } else if (atom_tmp.get_atom_name() == "C4' ") {
          tmpsugar_coors[1] = atom_tmp.get_coordinate();
        } else if (atom_tmp.get_atom_name() == "O4' ") {
          tmpsugar_coors[2] = atom_tmp.get_coordinate();
        } else if (atom_tmp.get_atom_name() == "C3' ") {
          tmpsugar_coors[3] = atom_tmp.get_coordinate();
        } else if (atom_tmp.get_atom_name() == "C2' ") {
          tmpsugar_coors[4] = atom_tmp.get_coordinate();
        } else if (atom_tmp.get_atom_name() == "C1' ") {
          tmpsugar_coors[5] = atom_tmp.get_coordinate();
        }
        tmp_atoms.push_back(atom_tmp);
      }
      tmp_sugar.set_conformation(tmpsugar_coors);
      while (1) {
        infile >> atom_tmp;
        if (atom_tmp.get_record_name() == "ENDMDL") {
          tmp_atoms.push_back(atom_tmp);
          break;
        }
        tmp_atoms.push_back(atom_tmp);
      }
      pinang::Transform t;
      pinang::find_transform(tmp_sugar, std_sugar, t);
      for (j = 0; j < tmp_atoms.size(); ++j) {
        atmp = tmp_atoms[j];
        if (atmp.get_record_name() != "ATOM  ") {
          outfile << atmp;
          continue;
        }
        vcoor = atmp.get_coordinate();
        vnew = t.apply(vcoor);
        atmp.set_coordinate(vnew);
        outfile << atmp;
      }
      tmp_atoms.clear();
    }
  }

  // -------------------- Calculating... --------------------

  infile.close();
  reffile.close();
  outfile.close();

  return 0;
}

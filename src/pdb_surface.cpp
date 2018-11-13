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
  int in_flag = 0;
  int out_flag = 0;

  string basefilename = "";
  string infilename = "some.pdb";
  string outfilename = "some.dat";

  double grid_padding = 8.0;
  double grid_size = 2.0;
  double tip_size = 4.0;

  while ((opt = getopt(argc, argv, "s:t:p:o:f:h")) != -1) {
    switch (opt) {
    case 's':
      grid_size = atof(optarg);
      break;
    case 't':
      tip_size = atof(optarg);
      break;
    case 'p':
      grid_padding = atof(optarg);
      break;
    case 'o':
      outfilename = optarg;
      out_flag = 1;
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

  double tip_size_sq = tip_size * tip_size;

  if (!in_flag)
    {
      cout << " ERROR: need parameter for option -f: " << "\n";
      print_usage(argv[0]);
    }

  if (!out_flag) {
    outfilename = basefilename + "_surface.dat";
  }
  ofstream out_file(outfilename.c_str());


  /////////////////////////////////////////////////////////////////////////////
  //                          Stupid local variables                         //
  /////////////////////////////////////////////////////////////////////////////
  pinang::PDB pdb1(infilename);
  pinang::Model m0 = pdb1.get_model(mod_index);
  double x, y, z;
  int i = 0;
  int j = 0;
  int k = 0;
  int n_chain = m0.get_size();
  int n_residue = 0;
  int n_atom = 0;
  pinang::Chain c_tmp;
  pinang::ChainType ct;
  pinang::Residue r_tmp;
  pinang::Atom a_tmp;
  pinang::Vec3d coor_tmp;
  double pdb_x_max=-10000, pdb_x_min=10000;
  double pdb_y_max=-10000, pdb_y_min=10000;
  double pdb_z_max=-10000, pdb_z_min=10000;

  /////////////////////////////////////////////////////////////////////////////
  //                         Flat structure of atoms                         //
  /////////////////////////////////////////////////////////////////////////////
  vector<pinang::Atom> f_atoms;
  vector<int> surface_residue_flag;

  // get pdb_max, pdb_min /////////////////////////////////////////////////////
  for (i = 0; i < n_chain; ++i) {
    c_tmp = m0.get_chain(i);
    ct = c_tmp.get_chain_type();
    if (ct == pinang::water || ct == pinang::other || ct == pinang::none)
      continue;
    n_residue = c_tmp.get_size();
    for (j = 0; j < n_residue; ++j) {
      r_tmp = c_tmp.get_residue(j);
      n_atom = r_tmp.get_size();
      for (k = 0; k < n_atom; ++k) {
        a_tmp = r_tmp.get_atom(k);
        coor_tmp = a_tmp.get_coordinate();
        x = coor_tmp.x();
        y = coor_tmp.y();
        z = coor_tmp.z();
        if (x > pdb_x_max) pdb_x_max = x;
        if (y > pdb_y_max) pdb_y_max = y;
        if (z > pdb_z_max) pdb_z_max = z;
        if (x < pdb_x_min) pdb_x_min = x;
        if (y < pdb_y_min) pdb_y_min = y;
        if (z < pdb_z_min) pdb_z_min = z;
        f_atoms.push_back(a_tmp);
        surface_residue_flag.push_back(0);
      }
    }
  }
  // cout << pdb_x_max << "  " << pdb_x_min << "\n";
  // cout << pdb_y_max << "  " << pdb_y_min << "\n";
  // cout << pdb_z_max << "  " << pdb_z_min << "\n";
  int n_f_atoms = f_atoms.size();
  cout << "Number of atoms:" << n_f_atoms << "\n";

  /////////////////////////////////////////////////////////////////////////////
  //                                Make Grid                                //
  /////////////////////////////////////////////////////////////////////////////
  int grid_n_x = 0, grid_n_y = 0, grid_n_z = 0;
  vector<double> grid_real_x;
  vector<double> grid_real_y;
  vector<double> grid_real_z;

  double grid_x_max = pdb_x_max + grid_padding;
  double grid_y_max = pdb_y_max + grid_padding;
  double grid_z_max = pdb_z_max + grid_padding;
  double grid_x_min = pdb_x_min - grid_padding;
  double grid_y_min = pdb_y_min - grid_padding;
  double grid_z_min = pdb_z_min - grid_padding;

  x = grid_x_min;
  while (x <= grid_x_max) {
    grid_real_x.push_back(x);
    ++grid_n_x;
    x += grid_size;
  }
  y = grid_y_min;
  while (y <= grid_y_max) {
    grid_real_y.push_back(y);
    ++grid_n_y;
    y += grid_size;
  }
  z = grid_z_min;
  while (z <= grid_z_max) {
    grid_real_z.push_back(z);
    ++grid_n_z;
    z += grid_size;
  }
  int grid_resolution = grid_n_x * grid_n_y * grid_n_z;
  cout << "Number of grids = " << grid_n_x << " * "
       << grid_n_y << " * "
       << grid_n_z << " = "
       << grid_resolution
       << "\n";


  /////////////////////////////////////////////////////////////////////////////
  //                          Pick surrounding grids                         //
  /////////////////////////////////////////////////////////////////////////////
  cout << "============================================================" << "\n";
  cout << " Determining surrounding grids ... " << "\n";
  vector<int> surrounding_grid_flag(grid_resolution, 1);
  int i_grid = 0;
  double dist_sq = 0;
  double grid_x = 0;
  double grid_y = 0;
  double grid_z = 0;
  double dx = 0, dy = 0, dz = 0;
  for (int ix=0; ix < grid_n_x; ++ix) {
    grid_x = grid_real_x[ix];
    for (int iy=0; iy < grid_n_y; ++iy) {
      grid_y = grid_real_y[iy];
      for (int iz=0; iz < grid_n_z; ++iz) {
        grid_z = grid_real_z[iz];

        for (i = 0; i < n_f_atoms; ++i) {
          a_tmp = f_atoms[i];
          coor_tmp = a_tmp.get_coordinate();
          x = coor_tmp.x();
          y = coor_tmp.y();
          z = coor_tmp.z();
          dx = x - grid_x;
          if (dx < -tip_size || dx > tip_size) continue;
          dy = y - grid_y;
          if (dy < -tip_size || dy > tip_size) continue;
          dz = z - grid_z;
          if (dz < -tip_size || dz > tip_size) continue;
          dist_sq = dx * dx + dy * dy + dz * dz;
          if (dist_sq <= tip_size_sq) {
            surrounding_grid_flag[i_grid] = 0;
            break;
          }
        }

        ++i_grid;
      } // iz
    } // iy
    cout << ".";
  } // ix
  cout << " Done! " << endl;

  /////////////////////////////////////////////////////////////////////////////
  //                        Find out surface residues                        //
  /////////////////////////////////////////////////////////////////////////////
  cout << "============================================================" << "\n";
  cout << " Determining surface residues ... " << "\n";

  i_grid = 0;
  for (int ix=0; ix < grid_n_x; ++ix) {
    grid_x = grid_real_x[ix];
    for (int iy=0; iy < grid_n_y; ++iy) {
      grid_y = grid_real_y[iy];
      for (int iz=0; iz < grid_n_z; ++iz) {
        grid_z = grid_real_z[iz];

        if (!surrounding_grid_flag[i_grid++])
          continue;

        double dist_min_sq = 1.0e20;
        int i_surf_atom = -1;
        for (i = 0; i < n_f_atoms; ++i) {
          a_tmp = f_atoms[i];
          coor_tmp = a_tmp.get_coordinate();
          x = coor_tmp.x();
          y = coor_tmp.y();
          z = coor_tmp.z();
          dx = x - grid_x;
          dy = y - grid_y;
          dz = z - grid_z;
          dist_sq = dx * dx + dy * dy + dz * dz;
          if (dist_min_sq > dist_sq) {
            dist_min_sq = dist_sq;
            i_surf_atom = i;
          }
        }
        surface_residue_flag[i_surf_atom] = 1;
      } // iz
    } // iy
    cout << ".";
  } // ix
  cout << " Done! " << endl;

  cout << "============================================================" << "\n";
  cout << " Output ... " << endl;
  char chainID_tmp = '?';
  string resName_tmp = "Whatever";
  int resSerial_tmp = -10000;
  for (i = 0; i < n_f_atoms; ++i) {
    if (surface_residue_flag[i]) {
      if (chainID_tmp == f_atoms[i].get_chain_ID()
          && resName_tmp.compare(f_atoms[i].get_residue_name()) == 0
          && resSerial_tmp == f_atoms[i].get_residue_serial()) 
        continue;
      chainID_tmp = f_atoms[i].get_chain_ID();
      resName_tmp = f_atoms[i].get_residue_name();
      resSerial_tmp = f_atoms[i].get_residue_serial();
      // cout << chainID_tmp << "   " << resSerial_tmp << "  " << resName_tmp << endl;
      out_file << chainID_tmp << "   " << resSerial_tmp << "  " << resName_tmp << endl;
    }
  }
  cout << " Finish! " << endl;
  cout << "============================================================" << "\n";

  out_file.close();

  return 0;
}

void print_usage(char* s)
{
  cout << " Usage: "
       << s
       << "\n\t -f some.pdb\n\t"
       << " [-o (out_surface_resid_index.dat)]\n\t"
       << " [-t tip_size]\n\t"
       << " [-s grid_size]\n\t"
       << " [-p grid_padding]\n\t"
       << " [-h]"
       << "\n";
  exit(EXIT_SUCCESS);
}

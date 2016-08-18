/*!
  @file model.cpp
  @brief Define functions of class Model.

  Definitions of member or friend functions of class Model.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:43
  @copyright GNU Public License V3.0
*/

#include <iomanip>
#include "model.hpp"

namespace pinang {

Chain& Model::get_chain(unsigned int n)
{
  if (v_chains_.empty())
  {
    std::cout << " ~             PINANG :: model.hpp              ~ " << "\n";
    std::cerr << "ERROR: No Chains found in Model: "
              << model_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  if (n >= v_chains_.size())
  {
    std::cout << " ~             PINANG :: model.hpp              ~ " << "\n";
    std::cerr << "ERROR: Chain number out of range in Model: "
              << model_serial_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return v_chains_[n];
}

void Model::add_chain(Chain& c)
{
  c.self_check();
  v_chains_.push_back(c);
  ++n_chain_;
}


void Model::output_sequence(int n) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    v_chains_[i].output_sequence(n);
}

void Model::output_sequence_fasta(std::ostream & f_fasta, std::string s) const
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i)
    v_chains_[i].output_sequence_fasta(f_fasta, s);
}

// Model ===================================================================
Model::Model()
{
  model_serial_ = 0;
  v_chains_.clear();
  n_chain_ = 0;
}

void Model::reset()
{
  model_serial_ = 0;
  v_chains_.clear();
  n_chain_ = 0;
}

void Model::output_cg_crd(std::ostream& o)
{
  int i = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_cg_crd(o);
  }
}

void Model::output_cg_pdb(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_cg_pdb(o, n, m);
  }
  o << "ENDMDL" << std::endl;
}

void Model::output_top_mass(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 3 - 1;
    else
      n += v_chains_[i].get_size();
  }

  o << std::setw(8) << n << " !NATOM \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none) {
      continue;
    }
    v_chains_[i].output_top_mass(o, n, m);
  }
  o << std::endl;
}

void Model::output_top_bond(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none || ct == ion)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 3 - 2;
    else
      n += v_chains_[i].get_size() - 1;
  }

  o << std::setw(8) << n << " !NBOND: bonds \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_top_bond(o, n, m);
  }
  o << " \n" << std::endl;
}

void Model::output_top_angle(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none || ct == ion)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 4 - 5;
    else
      n += v_chains_[i].get_size() - 2;
  }

  o << std::setw(8) << n << " !NTHETA: angles \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_top_angle(o, n, m);
  }
  o << " \n" << std::endl;
}

void Model::output_top_dihedral(std::ostream& o)
{
  int i = 0;
  int n = 0;
  int m = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none || ct == ion)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 2 - 4;
    else
      n += v_chains_[i].get_size()-3;
  }

  o << std::setw(8) << n << " !NPHI: dihedrals \n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_top_dihedral(o, n, m);
  }
  o << "\n" << std::endl;
}

void Model::output_ffparm_bond(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 3 - 2;
    else
      n += v_chains_[i].get_size() - 1;
  }

  o << "[ bonds ]" << std::setw(8) << n << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj"
    << std::setw(17) << "r0" << std::setw(9) << "K_b" << "\n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_ffparm_bond(o, n);
  }
  o << " \n" << std::endl;
}

void Model::output_ffparm_angle(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 4 - 5;
    else
      n += v_chains_[i].get_size()-2;
  }

  o << "[ angles ]" << std::setw(8) << n << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj" << std::setw(9) << "pk"
    << std::setw(13) << "theta_0" << std::setw(9) << "K_a" << "\n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_ffparm_angle(o, n);
  }
  o << "\n" << std::endl;
}

void Model::output_ffparm_dihedral(std::ostream& o)
{
  int i = 0;
  int n = 0;

  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    else if (ct == DNA || ct == RNA || ct == na)
      n += v_chains_[i].get_size() * 2 - 4;
    else
      n += v_chains_[i].get_size()-3;
  }

  o << "[ dihedrals ]" << std::setw(8) << n << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj"
    << std::setw(9) << "pk" << std::setw(9) << "pl"
    << std::setw(13) << "phi_0" << std::setw(9) << "K_d_1"
    << std::setw(9) << "K_d_3" << "\n";

  n = 0;
  for (i = 0; i < n_chain_; ++i) {
    ChainType ct = v_chains_[i].get_chain_type();
    if (ct == water || ct == other || ct == none)
      continue;
    v_chains_[i].output_ffparm_dihedral(o, n);
  }
  o << "\n" << std::endl;
}

void Model::output_ffparm_nonbonded(std::ostream& o)
{
  int i, j, k;
  Atom atmp1, atmp2, atmp3, atmp4;
  Residue rtmp1, rtmp2;
  ChainType ct_tmp;
  int pg_size, dg_size;
  double cg_dist = 0;
  double aa_dist_min = 0;
  std::string groove_info;

  std::vector<Atom> tmp_cg_pro_group;
  std::vector<Atom> tmp_cg_dna_group;
  std::vector<Residue> tmp_residue_pro_group;
  std::vector<Residue> tmp_residue_dna_group;

  int tmp_resid_serial = 0;
  int tmp_chain_serial = -1;
  for (i = 0; i < n_chain_; ++i) {
    ct_tmp = v_chains_[i].get_chain_type();
    if (ct_tmp == none || ct_tmp == water || ct_tmp == other)
      continue;
    ++tmp_chain_serial;
    int m_chain_size = v_chains_[i].get_size();
    if (ct_tmp == protein) {
      for (j = 0; j < m_chain_size; ++j) {
        rtmp1 = v_chains_[i].get_residue(j);
        atmp1 = rtmp1.get_cg_C_alpha();
        ++tmp_resid_serial;
        atmp1.set_residue_serial(tmp_resid_serial);
        atmp1.set_chain_ID(tmp_chain_serial + 97);
        tmp_cg_pro_group.push_back(atmp1);
        tmp_residue_pro_group.push_back(rtmp1);
      }
    } else if (ct_tmp == DNA) {
      Residue rtmp_P, rtmp_S, rtmp_B, rtmp_P_1;
      for (j = 0; j < m_chain_size; ++j) {
        // Add CG Phosphate;
        rtmp1 = v_chains_[i].get_residue(j);
        if (j != 0) {
          atmp1 = rtmp1.get_cg_P();
          ++tmp_resid_serial;
          atmp1.set_residue_serial(tmp_resid_serial);
          atmp1.set_chain_ID(tmp_chain_serial + 97);
          tmp_cg_dna_group.push_back(atmp1);
        }
        // Add CG Sugar;
        atmp1 = rtmp1.get_cg_S();
        ++tmp_resid_serial;
        atmp1.set_residue_serial(tmp_resid_serial);
        atmp1.set_chain_ID(tmp_chain_serial + 97);
        tmp_cg_dna_group.push_back(atmp1);

        // Add CG Base;
        atmp1 = rtmp1.get_cg_B();
        ++tmp_resid_serial;
        atmp1.set_residue_serial(tmp_resid_serial);
        atmp1.set_chain_ID(tmp_chain_serial + 97);
        tmp_cg_dna_group.push_back(atmp1);

        rtmp_P = rtmp_P_1;
        rtmp_P_1.reset();
        rtmp_P.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_P_1.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_S.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_B.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_P.set_chain_ID(rtmp1.get_chain_ID());
        rtmp_P_1.set_chain_ID(rtmp1.get_chain_ID());
        rtmp_S.set_chain_ID(rtmp1.get_chain_ID());
        rtmp_B.set_chain_ID(rtmp1.get_chain_ID());

        for (k = 0; k < rtmp1.get_size(); ++k) {
          atmp1 = rtmp1.get_atom(k);
          std::string aname = atmp1.get_atom_name();
          if (aname == "O3' ") {
            rtmp_P_1.add_atom(atmp1);
          } else if (aname == "P   " || aname == "OP1 " || aname == "OP2 " ||
                     aname == "O5' " || aname == "O1P " || aname == "O2P ") {
            rtmp_P.add_atom(atmp1);
          } else {
            std::string::size_type np;
            np = aname.find("'");
            if (np == std::string::npos) {
              rtmp_B.add_atom(atmp1);
            } else {
              rtmp_S.add_atom(atmp1);
            }
          }
        }
        if (j != 0)
          tmp_residue_dna_group.push_back(rtmp_P);
        tmp_residue_dna_group.push_back(rtmp_S);
        tmp_residue_dna_group.push_back(rtmp_B);
        rtmp_P.reset();
        rtmp_S.reset();
        rtmp_B.reset();
      }
    } else if (ct_tmp == RNA || ct_tmp == na) {
      for (j = 0; j < m_chain_size; ++j) {
        if (j != 0) {
          tmp_resid_serial += 3;
        } else {
          tmp_resid_serial += 2;
        }
      }
    } else if (ct_tmp == ion) {
      ++tmp_resid_serial;
    }
  }

  // Check vector sizes...
  if (tmp_cg_pro_group.size() != tmp_residue_pro_group.size()) {
    std::cout << " ~      PINANG::model_statistics.cpp     ~ \n";
    std::cout << "ERROR in getting protein atom group and residue group. \n";
    exit(EXIT_SUCCESS);
  } else {
    pg_size = int(tmp_cg_pro_group.size());
  }
  if (tmp_cg_dna_group.size() != tmp_residue_dna_group.size()) {
    std::cout << " ~      PINANG::model_statistics.cpp     ~ \n";
    std::cout << "ERROR in getting DNA atom group and residue group. \n";
    exit(EXIT_SUCCESS);
  } else {
    dg_size = int(tmp_cg_dna_group.size());
  }

  // Computing protein-protein native contacts...
  std::vector<int> pro_contact_part_1_atom_serial;
  std::vector<int> pro_contact_part_2_atom_serial;
  std::vector<int> pro_contact_part_1_chain_ID;
  std::vector<int> pro_contact_part_2_chain_ID;
  std::vector<double> pro_contact_cg_distance;
  for (i = 0; i < pg_size - 1; ++i) {
    atmp1 = tmp_cg_pro_group[i];
    rtmp1 = tmp_residue_pro_group[i];
    for (j = i + 1; j < pg_size; ++j) {
      atmp2 = tmp_cg_pro_group[j];
      rtmp2 = tmp_residue_pro_group[j];
      if (atmp2.get_residue_serial() <= atmp1.get_residue_serial() + 3 && atmp1.get_chain_ID() == atmp2.get_chain_ID())
        continue;
      aa_dist_min = residue_min_distance(rtmp1, rtmp2);
      if (aa_dist_min < g_pro_pro_aa_cutoff && aa_dist_min > 0) {
        cg_dist = atom_distance(atmp1, atmp2);
        pro_contact_part_1_atom_serial.push_back(atmp1.get_residue_serial());
        pro_contact_part_2_atom_serial.push_back(atmp2.get_residue_serial());
        pro_contact_part_1_chain_ID.push_back(atmp1.get_chain_ID());
        pro_contact_part_2_chain_ID.push_back(atmp2.get_chain_ID());
        pro_contact_cg_distance.push_back(cg_dist);
      }
    }
  }
  o << "[ native ]" << std::setw(8) << pro_contact_cg_distance.size() << "\n";
  o << "# " << std::setw(6) << "pi" << std::setw(9) << "pj"
    << std::setw(5) << "ci" << std::setw(5) << "cj"
    << std::setw(17) << "sigma" << std::setw(13) << "eps" << "\n";
  for (i = 0; i < pro_contact_cg_distance.size(); ++i) {
    o << std::setw(8) << pro_contact_part_1_atom_serial[i] << " "
      << std::setw(8) << pro_contact_part_2_atom_serial[i] << " "
      << std::setw(4) << char(pro_contact_part_1_chain_ID[i]) << " "
      << std::setw(4) << char(pro_contact_part_2_chain_ID[i]) << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
      << std::setw(16) << pro_contact_cg_distance[i] << " "
      << std::setprecision(5) << std::setw(12) << k_K_native << " " << "\n";
  }
  o << std::endl;

  // Computing protein-DNA contacts...
  pinang::Vec3d tmp_c_CA, tmp_c_B0, tmp_c_S0, tmp_c_B5, tmp_c_B3;
  pinang::Vec3d tmp_c_CA_N, tmp_c_CA_C;
  double tmp_angle_0, tmp_angle_53, tmp_angle_NC;
  int calpha_term_N = 0, calpha_term_C = 0;
  std::vector<int> pro_DNA_contact_pro_atom_serial;
  std::vector<int> pro_DNA_contact_DNA_atom_serial;
  std::vector<double> pro_DNA_contact_cg_distance;
  std::vector<double> pro_DNA_contact_cg_angle_NC;
  std::vector<double> pro_DNA_contact_cg_angle_0;
  std::vector<double> pro_DNA_contact_cg_angle_53;
  std::vector<std::string> pro_DNA_contact_cg_groove_info;
  for (i = 0; i < pg_size; ++i) {
    atmp1 = tmp_cg_pro_group[i];
    rtmp1 = tmp_residue_pro_group[i];
    tmp_c_CA = atmp1.get_coordinate();  // Coor of CA
    k = i - 1;
    if (k < 0 || tmp_cg_pro_group[k].get_chain_ID() != atmp1.get_chain_ID()) {
      tmp_c_CA_N = tmp_c_CA;
      calpha_term_N = 1;
    } else {
      calpha_term_N = 0;
      tmp_c_CA_N = tmp_cg_pro_group[k].get_coordinate();  // Coor of N' CA
    }
    k = i + 1;
    if (k >= tmp_cg_pro_group.size() || tmp_cg_pro_group[k].get_chain_ID() != atmp1.get_chain_ID()) {
      tmp_c_CA_C = tmp_c_CA;
      calpha_term_C = 1;
    } else {
      calpha_term_C = 0;
      tmp_c_CA_C = tmp_cg_pro_group[k].get_coordinate();  // Coor of C' CA
    }
    if (calpha_term_N * calpha_term_C > 0) {
      std::cout << " Single residue Chain!!! WTF!!! \n";
      exit(EXIT_SUCCESS);
    }
 
    for (j = 0; j < dg_size; ++j) {
      atmp2 = tmp_cg_dna_group[j];
      rtmp2 = tmp_residue_dna_group[j];
      if (atmp2.get_atom_name() != "DB  ")
        continue;
      aa_dist_min = residue_min_distance(rtmp1, rtmp2, atmp3, atmp4);
      if (aa_dist_min < g_pro_DNA_aa_cutoff && aa_dist_min > 0) {
        // ---------- 5' Base ----------
        k = j - 3;
        if (k < 0 || tmp_cg_dna_group[k].get_chain_ID() != tmp_cg_dna_group[j].get_chain_ID()) {
          continue;
        } else {
          if (tmp_cg_dna_group[k].get_atom_name() != "DB  ") {
            std::cout << " Wrong interaction pair for 5'DB!!! \n";
            exit(EXIT_SUCCESS);
          }
          tmp_c_B5 = tmp_cg_dna_group[k].get_coordinate();  // Coor of 5' B
        }
        // ---------- 3' Base ----------
        k = j + 3;
        if (k > tmp_cg_dna_group.size() || tmp_cg_dna_group[k].get_chain_ID() != tmp_cg_dna_group[j].get_chain_ID()) {
          continue;
        } else {
          if (tmp_cg_dna_group[k].get_atom_name() != "DB  ") {
            std::cout << " Wrong interaction pair for 3'DB!!! \n";
            exit(EXIT_SUCCESS);
          }
          tmp_c_B3 = tmp_cg_dna_group[k].get_coordinate();  // Coor of 5' B
        }
        tmp_c_B0 = atmp2.get_coordinate();  // Coor of B
        // ---------- Sugar -- Base -- CA angle ----------
        if (tmp_cg_dna_group[j - 1].get_atom_name() != "DS  ") {
          std::cout << " Wrong interaction pair for DS in DS-DB-CA!!! \n";
          exit(EXIT_SUCCESS);
        }
        tmp_c_S0 = tmp_cg_dna_group[j-1].get_coordinate();  // Coor of S connected to B
        // ---------------------------------------- CALCULATIONS! ----------------------------------------
        cg_dist = atom_distance(atmp1, atmp2);
        groove_info = get_DNA_atom_position_info(atmp4.get_residue_name(), atmp4.get_atom_name());
        tmp_angle_53 = vec_angle_deg(tmp_c_B3 - tmp_c_B5, tmp_c_CA - tmp_c_B0);
        tmp_angle_0 = vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_c_CA - tmp_c_B0);
        tmp_angle_NC = vec_angle_deg(tmp_c_CA_C - tmp_c_CA_N, tmp_c_B0 - tmp_c_CA);

        pro_DNA_contact_pro_atom_serial.push_back(atmp1.get_residue_serial());
        pro_DNA_contact_DNA_atom_serial.push_back(atmp2.get_residue_serial());
        pro_DNA_contact_cg_distance.push_back(cg_dist);
        pro_DNA_contact_cg_angle_NC.push_back(tmp_angle_NC);
        pro_DNA_contact_cg_angle_0.push_back(tmp_angle_0);
        pro_DNA_contact_cg_angle_53.push_back(tmp_angle_53);
        pro_DNA_contact_cg_groove_info.push_back(groove_info);
      }
    }
  }
  o << "[ protein-DNA seq-specific ]" << std::setw(8) << pro_DNA_contact_cg_distance.size() << "\n";
  o << "# " << std::setw(6) << "pro_i" << std::setw(9) << "dna_j"
    << std::setw(12) << "r_0" << std::setw(9) << "angle_NC"
    // << std::setw(9) << "angle_53" << std::setw(10) << "groove"
    << std::setw(9) << "angle_0" 
    << std::setw(9) << "angle_53" 
    << std::setw(9) << "sigma" << std::setw(9) << "phi" << "\n";
  for (i = 0; i < pro_DNA_contact_cg_distance.size(); ++i) {
    o << std::setw(8) << pro_DNA_contact_pro_atom_serial[i] << " "
      << std::setw(8) << pro_DNA_contact_DNA_atom_serial[i] << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
      << std::setw(11) << pro_DNA_contact_cg_distance[i] << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
      << std::setw(8) << pro_DNA_contact_cg_angle_NC[i] << " "
      << std::setw(8) << pro_DNA_contact_cg_angle_0[i] << " "
      << std::setw(8) << pro_DNA_contact_cg_angle_53[i] << " "
      // << std::setw(9) << pro_DNA_contact_cg_groove_info[i] << " "
      << std::setprecision(3) << std::setw(8)
      << 1.0 << " " << std::setw(8) << 10.0 << " " << "\n";
  }
  o << std::endl;
}

std::ostream& operator<<(std::ostream& o, Model& m)
{
  o << "MODEL "
    << std::setw(8) << m.model_serial_
    << "\n";
  int i = 0;
  int s = m.n_chain_;
  for (i = 0; i < s; ++i) {
    o << m.v_chains_[i];
  }
  o << "ENDMDL" << "\n";
  return o;
}

}  // pinang

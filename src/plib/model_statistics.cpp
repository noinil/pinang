/*!
  @file model_statistics.cpp
  @brief Statistics of structural properties.

  Statistics of structural properties such as distances between contacts of CG
  particles, angles formed between CG particles...

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-06-15 17:14
  @copyright GNU Public License V3.0
*/

#include <iomanip>
#include "model.hpp"

namespace pinang {

void Model::output_statistics_pro_DNA_contact_pairs(std::ostream& o)
{
  int i, j;
  Atom atmp1, atmp2, atmp3, atmp4;
  Residue rtmp1, rtmp2;
  ChainType ct_tmp;
  int pg_size, dg_size;
  double cg_dist = 0;
  double aa_dist_min = 0;
  const double aa_dist_cutoff = 6.5;

  std::vector<Atom> tmp_cg_pro_group;
  std::vector<Atom> tmp_cg_dna_group;
  std::vector<Residue> tmp_residue_pro_group;
  std::vector<Residue> tmp_residue_dna_group;

  int tmp_resid_serial = 0;
  for (i = 0; i < n_chain_; ++i) {
    ct_tmp = v_chains_[i].get_chain_type();
    if (ct_tmp == none || ct_tmp == water || ct_tmp == ion || ct_tmp == other || ct_tmp == RNA || ct_tmp == na)
      continue;
    int m_chain_size = v_chains_[i].get_size();
    if (ct_tmp == protein) {
      for (j = 0; j < m_chain_size; ++j) {
        rtmp1 = v_chains_[i].get_residue(j);
        atmp1 = rtmp1.get_cg_C_alpha();
        ++tmp_resid_serial;
        atmp1.set_residue_serial(tmp_resid_serial);
        tmp_cg_pro_group.push_back(atmp1);
        tmp_residue_pro_group.push_back(rtmp1);
      }
    } else if (ct_tmp == DNA) {
      Residue rtmp_P, rtmp_S, rtmp_B, rtmp_P_1;
      for (j = 0; j < m_chain_size; ++j) {
        // Add CG Phosphate;
        rtmp1 = v_chains_[i].get_residue(j);
        atmp1 = rtmp1.get_cg_P();
        if (j != 0) {
          ++tmp_resid_serial;
          atmp1.set_residue_serial(tmp_resid_serial);
          tmp_cg_dna_group.push_back(atmp1);
        }
        // Add CG Sugar;
        atmp1 = rtmp1.get_cg_S();
        ++tmp_resid_serial;
        atmp1.set_residue_serial(tmp_resid_serial);
        tmp_cg_dna_group.push_back(atmp1);

        // Add CG Base;
        atmp1 = rtmp1.get_cg_B();
        ++tmp_resid_serial;
        atmp1.set_residue_serial(tmp_resid_serial);
        tmp_cg_dna_group.push_back(atmp1);

        rtmp_P = rtmp_P_1;
        rtmp_P.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_S.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_B.set_residue_serial(rtmp1.get_residue_serial());
        rtmp_P_1.reset();
        rtmp_P_1.set_residue_serial(rtmp1.get_residue_serial() + 1);

        for (int k = 0; k < rtmp1.get_size(); ++k) {
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
    }
  }

  // Check vector sizes...
  if (tmp_cg_pro_group.size() != tmp_residue_pro_group.size()) {
    std::cout << " ~      PINANG::model_statistics.cpp     ~ \n";
    std::cout << "ERROR in getting protein atom group and residue group. \n";
    return;
  } else {
    pg_size = int(tmp_cg_pro_group.size());
  }
  if (tmp_cg_dna_group.size() != tmp_residue_dna_group.size()) {
    std::cout << " ~      PINANG::model_statistics.cpp     ~ \n";
    std::cout << "ERROR in getting DNA atom group and residue group. \n";
    return;
  } else {
    dg_size = int(tmp_cg_dna_group.size());
  }

  for (i = 0; i < pg_size; ++i) {
    atmp1 = tmp_cg_pro_group[i];
    rtmp1 = tmp_residue_dna_group[i];
    for (j = 0; j < dg_size; ++j) {
      atmp2 = tmp_cg_dna_group[j];
      rtmp2 = tmp_residue_dna_group[j];
      aa_dist_min = residue_min_distance(rtmp1, rtmp2, atmp3, atmp4);
      if (aa_dist_min < aa_dist_cutoff) {
        cg_dist = atom_distance(atmp1, atmp2);
        o << "aa : " << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
          << std::setw(9) << aa_dist_min << " : " << atmp3.get_chain_ID() << "|"
          << atmp3.get_residue_name() << " " << atmp3.get_residue_serial() << ":"
          << atmp3.get_atom_name() << "  --  " << atmp4.get_chain_ID() << "|"
          << atmp4.get_residue_name() << " " << atmp4.get_residue_serial() << ":"
          << atmp4.get_atom_name() << " \n";
        o << "cg : " << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
          << std::setw(9) << cg_dist << " : " << atmp1.get_chain_ID() << "|"
          << atmp1.get_residue_name() << " " << atmp1.get_residue_serial() << ":"
          << atmp1.get_atom_name() << "  --  " << atmp2.get_chain_ID() << "|"
          << atmp2.get_residue_name() << " " << atmp2.get_residue_serial() << ":"
          << atmp2.get_atom_name() << " \n";
      }
    }
  }

}

}

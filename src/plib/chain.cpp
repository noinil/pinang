/*!
  @file chain.cpp
  @brief Define functions of class Chain.

  Definitions of member or friend functions of class Chain.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-24 15:39
  @copyright GNU Public License V3.0
*/


#include <iomanip>
#include "chain.hpp"

namespace pinang {

Residue& Chain::get_residue(int n)
{
  if (v_residues_.empty())
  {
    std::cout << " ~             PINANG :: chain.hpp              ~ " << "\n";
    std::cerr << "ERROR: No Residues found in Chain: "
              << chain_ID_ << "\n";
    exit(EXIT_SUCCESS);
  }
  if (n >= int(v_residues_.size()))
  {
    std::cout << " ~             PINANG :: chain.hpp              ~ " << "\n";
    std::cerr << "ERROR: Residue index out of range in Chain: "
              << chain_ID_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return v_residues_[n];
}

int Chain::add_residue(const Residue& r)
{
  r.self_check();
  v_residues_.push_back(r);
  n_residue_++;
  return 0;
}

void Chain::self_check()
{
  for (Residue& r : v_residues_) {
    if (r.get_chain_ID() != chain_ID_ || r.get_chain_type() != chain_type_) {
      std::cout << " ~             PINANG :: chain.hpp              ~ " << "\n";
      std::cerr << "ERROR: Inconsistent chain ID or chain type in Chain "
                << chain_ID_ << "\n";
      exit(EXIT_SUCCESS);
    }
  }
  if (chain_type_ == protein) {
    v_residues_[0].set_terminus_flag(-1);
    v_residues_[n_residue_ - 1].set_terminus_flag(1);
  }
  if (chain_type_ == DNA) {
    v_residues_[0].set_terminus_flag(5);
    v_residues_[n_residue_ - 1].set_terminus_flag(3);
  }
}

void Chain::output_sequence(int n) const
{
  if (chain_type_ == water)
    return;

  int j = 0;

  std::cout << " - Chain " << chain_ID_
            << " : " << n_residue_ << " residues."
            << "\n";

  if (n == 1)
  {
    std::cout << " ";
    for (const Residue& r : v_residues_) {
      std::string s_tmp =  r.get_short_name();
      std::cout << std::setw(1) << s_tmp[1];
      j++;
      if (j%10 == 5)
        std::cout << " ";
      if (j%10 == 0)
      {
        std::cout << "  " << std::setw(4) << r.get_residue_serial() << "\n";
        std::cout << " ";
      }
    }
    std::cout << "\n";
  } else if (n == 3) {
    for (const Residue& r : v_residues_) {
      std::cout << std::setw(4) << r.get_residue_name();
      j++;
      if (j%10 == 0)
      {
        std::cout << "  " << std::setw(4) << r.get_residue_serial() << "\n";
      }
    }
    std::cout << "\n";
  }
}

void Chain::output_sequence_fasta(std::ostream & f_fasta, std::string s0) const
{
  if (chain_type_ == water)
    return;
  if (n_residue_ <= 3)
    return;

  f_fasta << ">" << s0 << "_chain_"
          << chain_ID_ << "_type_"
          << chain_type_
          << "\n";

  for (const Residue& r : v_residues_) {
    std::string s_tmp =  r.get_short_name();
    f_fasta << std::setw(1) << s_tmp[1];
  }
  f_fasta << "\n";
}

void Chain::output_cg_pos(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  o << " - Chain " << chain_ID_;
  o << " : " ;
  switch (chain_type_) {
    case protein:
      o << "protein";
      break;
    case DNA:
      o << "DNA";
      break;
    case RNA:
      o << "RNA";
      break;
    case ion:
      o << "ion";
      break;
    case na:
      o << "na";
      break;
    default:
      o << "unknown";
  }
  o << " : " << n_residue_
    << "\n";

  int i = 0;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (Residue& r : v_residues_) {
      Vec3d coor_CA;
      coor_CA = r.get_cg_C_alpha().get_coordinate();
      o << std::setw(8) << ++n
        << std::setw(5) << r.get_residue_name()
        << std::setw(5) << r.get_residue_serial() << "   "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
        << std::setw(12) << coor_CA.x() << " "
        << std::setw(12) << coor_CA.y() << " "
        << std::setw(12) << coor_CA.z() << " "
        << "\n";
    }
  } else {
    for (Residue& r : v_residues_) {
      Vec3d coor_P;
      Vec3d coor_S;
      Vec3d coor_B;
      if (r.get_terminus_flag() != 5) {
        coor_P = r.get_cg_P().get_coordinate();
        o << std::setw(8) << ++n
          << std::setw(5) << "P"
          << std::setw(5) << r.get_residue_serial() << "   "
          << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
          << std::setw(12) << coor_P.x() << " "
          << std::setw(12) << coor_P.y() << " "
          << std::setw(12) << coor_P.z() << " "
          << "\n";
      }
      coor_S = r.get_cg_S().get_coordinate();
      o << std::setw(8) << ++n
        << std::setw(5) << "S"
        << std::setw(5) << r.get_residue_serial() << "   "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
        << std::setw(12) << coor_S.x() << " "
        << std::setw(12) << coor_S.y() << " "
        << std::setw(12) << coor_S.z() << " "
        << "\n";
      coor_B = r.get_cg_B().get_coordinate();
      o << std::setw(8) << ++n
        << std::setw(5) << r.get_short_name()
        << std::setw(5) << r.get_residue_serial() << "   "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(3)
        << std::setw(12) << coor_B.x() << " "
        << std::setw(12) << coor_B.y() << " "
        << std::setw(12) << coor_B.z() << " "
        << "\n";
    }
  }
  o << "\n";
}

void Chain::output_top_mass(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
  {
    return;
  }

  int i = 0;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (const Residue& r : v_residues_) {
      o << std::setw(11) << ++n << " "
        << std::setw(8) << r.get_residue_serial() << " "
        << std::setw(9) << r.get_residue_name() << " "
        << std::setw(9) << "CA" << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
        << std::setw(16) << r.get_residue_mass() << " "
        << std::setw(12)
        << r.get_residue_charge()
        << "\n";
    }
  } else {
    for (const Residue& r : v_residues_) {
      if (r.get_terminus_flag() != 5) {
        o << std::setw(11) << ++n << " "
          << std::setw(8) << r.get_residue_serial() << " "
          << std::setw(9) << r.get_residue_name() << " "
          << std::setw(9) << "P" << " "
          << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
          << std::setw(16) << 94.93 << " "
          << std::setw(12) << -1.0
          << "\n";
      }
      o << std::setw(11) << ++n << " "
        << std::setw(8) << r.get_residue_serial() << " "
        << std::setw(9) << r.get_residue_name() << " "
        << std::setw(9) << "S" << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
        << std::setw(16) << 99.11 << " "
        << std::setw(12) << 0.0
        << "\n";
      o << std::setw(11) << ++n << " "
        << std::setw(8) << r.get_residue_serial() << " "
        << std::setw(9) << r.get_residue_name() << " "
        << std::setw(9) << "B" << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
        << std::setw(16) << r.get_residue_mass() - 94.93 - 99.11 << " "
        << std::setw(12) << 0.0
        << "\n";
    }
  }
}

void Chain::output_top_bond(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
  {
    return;
  }

  int i = 0;
  double d = 0;
  double d_ps = 0;
  double d_sb = 0;
  double d_sp = 0;

  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_ - 1; i++) {
      d = residue_ca_distance(v_residues_[i], v_residues_[i+1]);
      n++;
      o << std::setw(8) << n << " "
        << std::setw(8) << n + 1 << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6)
        << std::setw(16) << d << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_bond << " "
        << "\n";
    }
    n++;
  } else {
    d_sb = atom_distance(v_residues_[0].get_cg_S(), v_residues_[0].get_cg_B());
    o << std::setw(8) << n+1 << " "
      << std::setw(8) << n+2 << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
      << std::setw(16) << d_sb << " "
      << std::setprecision(1) << std::setw(8) << k_K_bond << " "
      << "\n";
    d_sp = atom_distance(v_residues_[1].get_cg_P(), v_residues_[0].get_cg_S());
    o << std::setw(8) << n+1 << " "
      << std::setw(8) << n+3 << " "
      << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(6)
      << std::setw(16) << d_sp << " "
      << std::setprecision(1)
      << std::setw(8) << k_K_bond << " "
      << "\n";
    for (i = 1; i < n_residue_ - 1; i++) {
      n += 3;
      d_ps = atom_distance(v_residues_[i].get_cg_P(), v_residues_[i].get_cg_S());
      o << std::setw(8) << n << " "
        << std::setw(8) << n+1 << " " << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6) << std::setw(16) << d_ps << " "
        << std::setprecision(1) << std::setw(8) << k_K_bond << " "
        << "\n";
      d_sb = atom_distance(v_residues_[i].get_cg_S(), v_residues_[i].get_cg_B());
      o << std::setw(8) << n+1 << " "
        << std::setw(8) << n+2 << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6) << std::setw(16) << d_sb << " "
        << std::setprecision(1) << std::setw(8) << k_K_bond << " "
        << "\n";
      d_sp = atom_distance(v_residues_[i].get_cg_S(), v_residues_[i+1].get_cg_P());
      o << std::setw(8) << n+1 << " "
        << std::setw(8) << n+3 << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6) << std::setw(16) << d_sp << " "
        << std::setprecision(1) << std::setw(8) << k_K_bond << " "
        << "\n";
    }
    i = n_residue_ - 1;
    n += 3;
    d_ps = atom_distance(v_residues_[i].get_cg_P(), v_residues_[i].get_cg_S());
    o << std::setw(8) << n << " "
      << std::setw(8) << n+1 << " "
      << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(6) << std::setw(16) << d_ps << " "
      << std::setprecision(1) << std::setw(8) << k_K_bond << " "
      << "\n";
    d_sb = atom_distance(v_residues_[i].get_cg_S(), v_residues_[i].get_cg_B());
    o << std::setw(8) << n+1 << " "
      << std::setw(8) << n+2 << " "
      << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(6) << std::setw(16) << d_sb << " "
      << std::setprecision(1) << std::setw(8) << k_K_bond << " "
      << "\n";
    n += 2;
  }
}

void Chain::output_top_angle(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double a = 0;
  Vec3d v1, v2;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-2; i++) {
      v1 = v_residues_[i].get_cg_C_alpha().get_coordinate()
           - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      v2 = v_residues_[i+2].get_cg_C_alpha().get_coordinate()
           - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      a = vec_angle_deg (v1, v2);
      o << std::setw(8) << i+1+n << " "
        << std::setw(8) << i+2+n << " "
        << std::setw(8) << i+3+n << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6)
        << std::setw(12) << a << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_angle
        << "\n";
    }
    n += n_residue_;
  } else {
    // ---------- angle BSP ----------
    v1 = v_residues_[0].get_cg_B().get_coordinate()
         - v_residues_[0].get_cg_S().get_coordinate();
    v2 = v_residues_[1].get_cg_P().get_coordinate()
         - v_residues_[0].get_cg_S().get_coordinate();
    a = vec_angle_deg (v1, v2);
    o << std::setw(8) << n+2 << " "
      << std::setw(8) << n+1 << " "
      << std::setw(8) << n+3 << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
      << std::setw(12) << a << " "
      << std::setprecision(1)
      << std::setw(8) << k_K_angle << "\n";
    // ---------- angle SPS ----------
    v1 = v_residues_[0].get_cg_S().get_coordinate()
         - v_residues_[1].get_cg_P().get_coordinate();
    v2 = v_residues_[1].get_cg_S().get_coordinate()
         - v_residues_[1].get_cg_P().get_coordinate();
    a = vec_angle_deg (v1, v2);
    o << std::setw(8) << n+1 << " "
      << std::setw(8) << n+3 << " "
      << std::setw(8) << n+4 << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
      << std::setw(12) << a << " "
      << std::setprecision(1)
      << std::setw(8) << k_K_angle << "\n";

    // -------------------- loop --------------------
    for (i = 1; i < n_residue_-1; i++) {
      n += 3;
      // ---------- angle PSB ----------
      v1 = v_residues_[i].get_cg_P().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i].get_cg_B().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      a = vec_angle_deg (v1, v2);
      o << std::setw(8) << n << " "
        << std::setw(8) << n+1 << " "
        << std::setw(8) << n+2 << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
        << std::setw(12) << a << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_angle << "\n";
      // ---------- angle PSP ----------
      v1 = v_residues_[i].get_cg_P().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i+1].get_cg_P().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      a = vec_angle_deg (v1, v2);
      o << std::setw(8) << n << " "
        << std::setw(8) << n+1 << " "
        << std::setw(8) << n+3 << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
        << std::setw(12) << a << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_angle << "\n";
      // ---------- angle BSP ----------
      v1 = v_residues_[i].get_cg_B().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i+1].get_cg_P().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      a = vec_angle_deg (v1, v2);
      o << std::setw(8) << n+2 << " "
        << std::setw(8) << n+1 << " "
        << std::setw(8) << n+3 << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
        << std::setw(12) << a << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_angle << "\n";
      // ---------- angle SPS ----------
      v1 = v_residues_[i].get_cg_S().get_coordinate()
           - v_residues_[i+1].get_cg_P().get_coordinate();
      v2 = v_residues_[i+1].get_cg_S().get_coordinate()
           - v_residues_[i+1].get_cg_P().get_coordinate();
      a = vec_angle_deg (v1, v2);
      o << std::setw(8) << n+1 << " "
        << std::setw(8) << n+3 << " "
        << std::setw(8) << n+4 << " "
        << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
        << std::setw(12) << a << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_angle << "\n";
    }
    n += 3;
    i = n_residue_ - 1;
    // ---------- angle PSB ----------
    v1 = v_residues_[i].get_cg_P().get_coordinate()
         - v_residues_[i].get_cg_S().get_coordinate();
    v2 = v_residues_[i].get_cg_B().get_coordinate()
         - v_residues_[i].get_cg_S().get_coordinate();
    a = vec_angle_deg (v1, v2);
    o << std::setw(8) << n << " "
      << std::setw(8) << n+1 << " "
      << std::setw(8) << n+2 << " "
      << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
      << std::setw(12) << a << " "
      << std::setprecision(1)
      << std::setw(8) << k_K_angle << "\n";
    n += 2;
  }
}

void Chain::output_top_dihedral(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double d = 0;           // dihedral
  Vec3d v1, v2, v3, n1, n2;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-3; i++) {
      v1 = v_residues_[i].get_cg_C_alpha().get_coordinate()
           - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      v2 = v_residues_[i+2].get_cg_C_alpha().get_coordinate()
           - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      v3 = v_residues_[i+2].get_cg_C_alpha().get_coordinate()
           - v_residues_[i+3].get_cg_C_alpha().get_coordinate();
      n1 = v1 % v2;
      n2 = v2 % v3;
      d = vec_angle_deg (n1, n2);
      o << std::setw(8) << i+1+n << " "
        << std::setw(8) << i+2+n << " "
        << std::setw(8) << i+3+n << " "
        << std::setw(8) << i+4+n << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6)
        << std::setw(12) << d << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_dihedral_1 << " "
        << std::setw(8) << k_K_dihedral_3
        << "\n";
    }
    n += n_residue_;
  } else {
    if (n_residue_ <= 2) {
      n += 5;
      return;
    }
    v1 = v_residues_[0].get_cg_S().get_coordinate()
         - v_residues_[1].get_cg_P().get_coordinate();
    v2 = v_residues_[1].get_cg_S().get_coordinate()
         - v_residues_[1].get_cg_P().get_coordinate();
    v3 = v_residues_[1].get_cg_S().get_coordinate()
         - v_residues_[2].get_cg_P().get_coordinate();

    n1 = v1 % v2;
    n2 = v2 % v3;
    d = vec_angle_deg (n1, n2);
    o << std::setw(8) << n + 1 << " "
      << std::setw(8) << n + 3 << " "
      << std::setw(8) << n + 4 << " "
      << std::setw(8) << n + 6 << " "
      << std::setiosflags(std::ios_base::fixed)
      << std::setprecision(6)
      << std::setw(12) << d << " "
      << std::setprecision(1)
      << std::setw(8) << k_K_dihedral_1 << " "
      << std::setw(8) << k_K_dihedral_3
      << "\n";

    for (i = 1; i < n_residue_-1; i++) {
      n += 3;
      // ---------- PSPS ----------
      v1 = v_residues_[i].get_cg_P().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i+1].get_cg_P().get_coordinate()
           - v_residues_[i].get_cg_S().get_coordinate();
      v3 = v_residues_[i+1].get_cg_P().get_coordinate()
           - v_residues_[i+1].get_cg_S().get_coordinate();
      n1 = v1 % v2;
      n2 = v2 % v3;
      d = vec_angle_deg (n1, n2);
      o << std::setw(8) << n << " "
        << std::setw(8) << n + 1 << " "
        << std::setw(8) << n + 3 << " "
        << std::setw(8) << 4 + n << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6)
        << std::setw(12) << d << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_dihedral_1 << " "
        << std::setw(8) << k_K_dihedral_3
        << "\n";

      if (i == n_residue_ - 2)
        break;
      // ---------- SPSP ----------
      v1 = v_residues_[i].get_cg_S().get_coordinate()
           - v_residues_[i+1].get_cg_P().get_coordinate();
      v2 = v_residues_[i+1].get_cg_S().get_coordinate()
           - v_residues_[i+1].get_cg_P().get_coordinate();
      v3 = v_residues_[i+1].get_cg_S().get_coordinate()
           - v_residues_[i+2].get_cg_P().get_coordinate();
      n1 = v1 % v2;
      n2 = v2 % v3;
      d = vec_angle_deg (n1, n2);
      o << std::setw(8) << n + 1 << " "
        << std::setw(8) << n + 3 << " "
        << std::setw(8) << n + 4 << " "
        << std::setw(8) << n + 6 << " "
        << std::setiosflags(std::ios_base::fixed)
        << std::setprecision(6)
        << std::setw(12) << d << " "
        << std::setprecision(1)
        << std::setw(8) << k_K_dihedral_1 << " "
        << std::setw(8) << k_K_dihedral_3
        << "\n";
    }
    n += 5;
  }
}

void Chain::output_top_native(std::ostream& o)
{
  int i = 0, j = 0;
  double d = -1, f = -1;
  ChainType cti, ctj;
  for (i = 0; i < n_residue_-4; i++) {
    cti = v_residues_[i].get_chain_type();
    if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
      continue;
    for (j = i + 4; j < n_residue_; j++) {
      ctj = v_residues_[j].get_chain_type();
      if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
        continue;
      d = residue_min_distance(v_residues_[i], v_residues_[j]);
      if ( d < g_cutoff)
      {
        f = residue_ca_distance(v_residues_[i], v_residues_[j]);
        o << std::setw(8) << i+1 << " "
          << std::setw(8) << j+1 << " "
          << std::setiosflags(std::ios_base::fixed)
          << std::setprecision(6)
          << std::setw(16) << f << " "
          << std::setprecision(5)
          << std::setw(12) << k_K_native << " "
          << "\n";
      }
    }
  }
}

int Chain::get_native_contact_number()
{
  int i = 0, j = 0;
  double d = -1;
  int n = 0;
  ChainType cti, ctj;

  for (i = 0; i < n_residue_-4; i++) {
    cti = v_residues_[i].get_chain_type();
    if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
      continue;
    for (j = i + 4; j < n_residue_; j++) {
      ctj = v_residues_[j].get_chain_type();
      if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
        continue;
      d = residue_min_distance(v_residues_[i], v_residues_[j]);
      if ( d < g_cutoff)
        n++;
    }
  }
  return n;
}

// Chain -------------------------------------------------------------------
Chain::Chain()
{
  chain_ID_ = -1;
  chain_type_ = none;
  v_residues_.clear();
  n_residue_ = 0;
}

void Chain::reset()
{
  chain_ID_ = -1;
  chain_type_ = none;
  v_residues_.clear();
  n_residue_ = 0;
}

Chain operator+(const Chain& c1, const Chain& c2)
{
  int i = 0;
  Chain c0;
  c0.set_chain_ID(c2.chain_ID_);
  int s1 = c1.n_residue_;
  int s2 = c2.n_residue_;
  if (s1 > 0)
    for (i = 0; i < s1; i++) {
      c0.v_residues_.push_back(c1.v_residues_[i]);
      c0.n_residue_++;
    }
  if (s2 > 0)
    for (i = 0; i < s2; i++) {
      c0.v_residues_.push_back(c2.v_residues_[i]);
      c0.n_residue_++;
    }

  return c0;
}

std::ostream& operator<<(std::ostream& o, Chain& c)
{
  int i = 0;
  int s = c.n_residue_;
  for (i = 0; i < s; i++) {
    o << c.v_residues_[i];
  }
  o << "TER   " << "\n";
  return o;
}

}  // pinang

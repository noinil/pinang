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
    std::cout << " ~             PINANG :: chain.hpp            ~ " << "\n";
    std::cerr << "ERROR: No Residues found in Chain: " << chain_ID_ << "\n";
    exit(EXIT_SUCCESS);
  }
  if (n >= int(v_residues_.size()))
  {
    std::cout << " ~             PINANG :: chain.hpp        ~ " << "\n";
    std::cerr << "ERROR: Residue index out of range in Chain: " << chain_ID_ << "\n";
    exit(EXIT_SUCCESS);
  }
  return v_residues_[n];
}

int Chain::add_residue(const Residue& r)
{
  r.self_check();
  v_residues_.push_back(r);
  ++n_residue_;
  return 0;
}

double get_mass_from_atom_name_tmp(const std::string& s)
{
  char c = s[0];
  switch (c) {
    case 'C':
      return 12.011;
    case 'O':
      return 15.999;
    case 'N':
      return 14.001;
    case 'P':
      return 30.974;
    case 'H':
      return 1.008;
    default:
      return 0.000;
  }
}
void Chain::self_check()
{
  for (Residue& r : v_residues_) {
    if (r.get_chain_ID() != chain_ID_ || r.get_chain_type() != chain_type_) {
      std::cout << " ~             PINANG :: chain.hpp              ~ " << "\n";
      std::cerr << "ERROR: Inconsistent chain ID or type in Chain " << chain_ID_ << "\n";
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

    // Recalculate CG coordinates for DNA particles.
    Vec3d com_P(0.0, 0.0, 0.0);
    Vec3d com_S(0.0, 0.0, 0.0);
    Vec3d com_B(0.0, 0.0, 0.0);
    double mass_P = 0.0;
    double mass_S = 0.0;
    double mass_B = 0.0;
    double mass_tmp_atom = 0.0;
    Vec3d tmp_O3p(0.0, 0.0, 0.0);
    int O3p_flag = 0;
    for (int i = 0; i < n_residue_; ++i) {
      Residue& r = v_residues_[i];
      if (r.get_cg_S().get_atom_name() == "DS  " && r.get_cg_B().get_atom_name() == "DB  ")
        continue;

      if (O3p_flag == 1) {
        com_P += tmp_O3p * 15.999;
        mass_P += 15.999;
      }
      O3p_flag = 0;
      int natom = r.get_size();
      for (int j = 0; j < natom; ++j) {
        Atom ta = r.get_atom(j);
        std::string aname = ta.get_atom_name();
        mass_tmp_atom = get_mass_from_atom_name_tmp(aname);
        if (aname == "O3' ") {
          tmp_O3p = ta.get_coordinate();
          O3p_flag = 1;
        } else if (aname == "P   " || aname == "OP1 " || aname == "OP2 " || aname == "O5' ") {
          com_P += ta.get_coordinate() * mass_tmp_atom;
          mass_P += mass_tmp_atom;
        } else {
          std::string::size_type np;
          np = aname.find("'");
          if (np == std::string::npos) {
            com_B += ta.get_coordinate() * mass_tmp_atom;
            mass_B += mass_tmp_atom;
          } else {
            com_S += ta.get_coordinate() * mass_tmp_atom;
            mass_S += mass_tmp_atom;
          }
        }
      }
      if (i > 0) {
        if (mass_P > 94.5) {
          com_P /= mass_P;
          r.get_cg_P().set_coordinate(com_P);
          r.get_cg_P().set_atom_name("DP");
        } else {
          std::cout << " ~             PINANG :: chain.hpp            ~ " << "\n";
          std::cerr << "ERROR: missing DNA atoms in Chain: " << chain_ID_ << "\n";
          std::cerr << "      in Residue: " << i << " 's phosphate group. \n";
          exit(EXIT_SUCCESS);
        }
      }
      if (mass_S > 75.5) {
        com_S /= mass_S;
        r.get_cg_S().set_coordinate(com_S);
        r.get_cg_S().set_atom_name("DS");
      } else {
        std::cout << " ~             PINANG :: chain.hpp            ~ " << "\n";
        std::cerr << "ERROR: missing DNA atoms in Chain: " << chain_ID_ << "\n";
        std::cerr << "      in Residue: " << i << " 's sugar group. \n";
        exit(EXIT_SUCCESS);
      }
      if (mass_B > r.get_residue_mass() - 94.97 - 83.11 - 8.0) {
        com_B /= mass_B;
        r.get_cg_B().set_coordinate(com_B);
        r.get_cg_B().set_atom_name("DB");
      } else {
        std::cout << " ~             PINANG :: chain.hpp            ~ " << "\n";
        std::cerr << "ERROR: missing DNA atoms in Chain: " << chain_ID_ << "\n";
        std::cerr << "      in Residue: " << i << " 's base. \n";
        exit(EXIT_SUCCESS);
      }
      mass_P = 0.0;
      mass_S = 0.0;
      mass_B = 0.0;
      com_P *= 0.0;
      com_S *= 0.0;
      com_B *= 0.0;
    }
  }
}

void Chain::output_sequence(int n) const
{
  if (chain_type_ == water)
    return;
  int j = 0;
  std::cout << " - Chain " << chain_ID_ << " : " << n_residue_ << " residues.\n";

  if (n == 1)
  {
    std::cout << " ";
    for (const Residue& r : v_residues_) {
      std::string s_tmp =  r.get_short_name();
      std::cout << std::setw(1) << s_tmp[1];
      ++j;
      if (j%10 == 5)
        std::cout << " ";
      if (j%10 == 0)
      {
        std::cout << "  " << std::setw(4) << r.get_residue_serial() << "\n ";
      }
    }
    std::cout << "\n";
  } else if (n == 3) {
    for (const Residue& r : v_residues_) {
      std::cout << std::setw(4) << r.get_residue_name();
      ++j;
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

  f_fasta << ">" << s0 << "_chain_" << chain_ID_ << "_type_" << chain_type_ << "\n";

  for (const Residue& r : v_residues_) {
    std::string s_tmp =  r.get_short_name();
    f_fasta << std::setw(1) << s_tmp[1];
  }
  f_fasta << "\n";
}

void Chain::output_cg_crd(std::ostream& o, int& n, int& m)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (Residue& r : v_residues_) {
      Atom pseudo_ca = r.get_cg_C_alpha();
      pseudo_ca.set_atom_serial(++n);
      pseudo_ca.set_chain_ID(m + 97);
      o << pseudo_ca;
    }
  } else {
    for (Residue& r : v_residues_) {
      Vec3d coor_B;
      if (r.get_terminus_flag() != 5) {
        Atom pseudo_P = r.get_cg_P();
        pseudo_P.set_atom_serial(++n);
        // pseudo_P.set_atom_name("DP");
        pseudo_P.set_chain_ID(m + 97);
        o << pseudo_P;
      }
      Atom pseudo_S = r.get_cg_S();
      pseudo_S.set_atom_serial(++n);
      // pseudo_S.set_atom_name("DS");
      pseudo_S.set_chain_ID(m + 97);
      o << pseudo_S;

      Atom pseudo_B = r.get_cg_B();
      pseudo_B.set_atom_serial(++n);
      // pseudo_B.set_atom_name("DB");
      pseudo_B.set_chain_ID(m + 97);
      o << pseudo_B;
    }
  }
  ++m;
  o << "TER\n";
}

void output_top_mass_line(std::ostream& o, int i, char s, int r, std::string rn,
                          std::string an, std::string at, double c, double m)
{
  o << std::setw(8) << i << std::setw(3) << s  << " " << std::setw(6) << r
    << " " << std::setw(3) << rn
    << std::setw(4) << an << " " << std::setw(4) << at << " "
    << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
    << std::setw(12) << c << " " << std::setprecision(6)
    << std::setw(13) << m << "           " << "0 \n";
}
void Chain::output_top_mass(std::ostream& o, int& n, int& m)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  char cid = m + 97;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (const Residue& r : v_residues_) {
      output_top_mass_line(o, ++n, cid, r.get_residue_serial(), r.get_residue_name(),
                           "CA", "C", r.get_residue_charge(), r.get_residue_mass());
    }
  } else {
    for (const Residue& r : v_residues_) {
      if (r.get_terminus_flag() != 5) {
        output_top_mass_line(o, ++n, cid, r.get_residue_serial(), r.get_residue_name(),
                             "DP", "P", -0.6, 94.97);
      }
      output_top_mass_line(o, ++n, cid, r.get_residue_serial(), r.get_residue_name(),
                           "DS", "S", 0.0, 83.11);
      output_top_mass_line(o, ++n, cid, r.get_residue_serial(), r.get_residue_name(),
                           "DB", "B", 0.0, r.get_residue_mass() - 83.11 - 94.97);
    }
  }
  ++m;
}

void Chain::output_top_bond(std::ostream& o, int& n, int& m)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  std::vector<int> out_bond_list;

  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_ - 1; ++i) {
      out_bond_list.push_back(++n);
      out_bond_list.push_back(n + 1);
    }
    ++n;
  } else {
    out_bond_list.push_back(n + 1);
    out_bond_list.push_back(n + 2);
    out_bond_list.push_back(n + 1);
    out_bond_list.push_back(n + 3);
    for (i = 1; i < n_residue_ - 1; ++i) {
      n += 3;
      out_bond_list.push_back(n);
      out_bond_list.push_back(n + 1);
      out_bond_list.push_back(n + 1);
      out_bond_list.push_back(n + 2);
      out_bond_list.push_back(n + 1);
      out_bond_list.push_back(n + 3);
    }
    n += 3;
    out_bond_list.push_back(n);
    out_bond_list.push_back(n + 1);
    out_bond_list.push_back(n + 1);
    out_bond_list.push_back(n + 2);
    n += 2;
  }

  for (int j : out_bond_list) {
    o << std::setw(8) << j;
    if (++m == 8) {
      o << "\n";
      m = 0;
    }
  }
}

void Chain::output_top_angle(std::ostream& o, int& n, int& m)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  std::vector<int> out_angle_list;
  int i = 0;

  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-2; ++i) {
      int ni = i + n;
      out_angle_list.push_back(ni + 1);
      out_angle_list.push_back(ni + 2);
      out_angle_list.push_back(ni + 3);
    }
    n += n_residue_;
  } else {
    // ---------- angle BSP ----------
    out_angle_list.push_back(n + 2);
    out_angle_list.push_back(n + 1);
    out_angle_list.push_back(n + 3);
    // ---------- angle SPS ----------
    out_angle_list.push_back(n + 1);
    out_angle_list.push_back(n + 3);
    out_angle_list.push_back(n + 4);

    // -------------------- loop --------------------
    for (i = 1; i < n_residue_-1; ++i) {
      n += 3;
      // ---------- angle PSB ----------
      out_angle_list.push_back(n);
      out_angle_list.push_back(n + 1);
      out_angle_list.push_back(n + 2);
      // ---------- angle PSP ----------
      out_angle_list.push_back(n);
      out_angle_list.push_back(n + 1);
      out_angle_list.push_back(n + 3);
      // ---------- angle BSP ----------
      out_angle_list.push_back(n + 2);
      out_angle_list.push_back(n + 1);
      out_angle_list.push_back(n + 3);
      // ---------- angle SPS ----------
      out_angle_list.push_back(n + 1);
      out_angle_list.push_back(n + 3);
      out_angle_list.push_back(n + 4);
    }
    n += 3;
    // ---------- angle PSB ----------
    out_angle_list.push_back(n);
    out_angle_list.push_back(n + 1);
    out_angle_list.push_back(n + 2);
    n += 2;
  }
  for (int j : out_angle_list) {
    o << std::setw(8) << j;
    if (++m == 9) {
      o << "\n";
      m = 0;
    }
  }
}

void Chain::output_top_dihedral(std::ostream& o, int& n, int& m)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  std::vector<int> out_dih_list;

  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-3; ++i) {
      int ni = i + n;
      out_dih_list.push_back(ni + 1);
      out_dih_list.push_back(ni + 2);
      out_dih_list.push_back(ni + 3);
      out_dih_list.push_back(ni + 4);
    }
    n += n_residue_;
  } else {
    if (n_residue_ <= 2) {
      n += 5;
      return;
    }
    out_dih_list.push_back(n + 1);
    out_dih_list.push_back(n + 3);
    out_dih_list.push_back(n + 4);
    out_dih_list.push_back(n + 6);

    for (i = 1; i < n_residue_-1; ++i) {
      n += 3;
      // ---------- PSPS ----------
      out_dih_list.push_back(n);
      out_dih_list.push_back(n + 1);
      out_dih_list.push_back(n + 3);
      out_dih_list.push_back(n + 4);

      if (i == n_residue_ - 2)
        break;
      // ---------- SPSP ----------
      out_dih_list.push_back(n + 1);
      out_dih_list.push_back(n + 3);
      out_dih_list.push_back(n + 4);
      out_dih_list.push_back(n + 6);
    }
    n += 5;
  }
  for (int j : out_dih_list) {
    o << std::setw(8) << j;
    if (++m == 8) {
      o << "\n";
      m = 0;
    }
  }
}

int Chain::get_native_contact_number()
{
  int i = 0, j = 0;
  double d = -1;
  int n = 0;
  ChainType cti, ctj;

  for (i = 0; i < n_residue_-4; ++i) {
    cti = v_residues_[i].get_chain_type();
    if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
      continue;
    for (j = i + 4; j < n_residue_; ++j) {
      ctj = v_residues_[j].get_chain_type();
      if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
        continue;
      d = residue_min_distance(v_residues_[i], v_residues_[j]);
      if ( d < g_cutoff)
        ++n;
    }
  }
  return n;
}

void output_ff_bond_line(std::ostream& o, int i, int j, double d, double k)
{
  o << std::setw(8) << i << " " << std::setw(8) << j << " "
    << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
    << std::setw(16) << d << " " << std::setprecision(1)
    << std::setw(8) << k << " \n";
}
void Chain::output_ffparm_bond(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double d = 0;
  double d_ps = 0;
  double d_sb = 0;
  double d_sp = 0;

  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_ - 1; ++i) {
      d = residue_ca_distance(v_residues_[i], v_residues_[i+1]);
      ++n;
      output_ff_bond_line(o, n, n + 1, d, k_K_bond);
    }
    ++n;
  } else {
    d_sb = atom_distance(v_residues_[0].get_cg_S(), v_residues_[0].get_cg_B());
    output_ff_bond_line(o, n + 1, n + 2, d_sb, k_K_bond);
    d_sp = atom_distance(v_residues_[1].get_cg_P(), v_residues_[0].get_cg_S());
    output_ff_bond_line(o, n + 1, n + 3, d_sp, k_K_bond);
    for (i = 1; i < n_residue_ - 1; ++i) {
      n += 3;
      d_ps = atom_distance(v_residues_[i].get_cg_P(), v_residues_[i].get_cg_S());
      output_ff_bond_line(o, n, n + 1, d_ps, k_K_bond);
      d_sb = atom_distance(v_residues_[i].get_cg_S(), v_residues_[i].get_cg_B());
      output_ff_bond_line(o, n + 1, n + 2, d_sb, k_K_bond);
      d_sp = atom_distance(v_residues_[i].get_cg_S(), v_residues_[i+1].get_cg_P());
      output_ff_bond_line(o, n + 1, n + 3, d_sp, k_K_bond);
    }
    i = n_residue_ - 1;
    n += 3;
    d_ps = atom_distance(v_residues_[i].get_cg_P(), v_residues_[i].get_cg_S());
    output_ff_bond_line(o, n, n + 1, d_ps, k_K_bond);
    d_sb = atom_distance(v_residues_[i].get_cg_S(), v_residues_[i].get_cg_B());
    output_ff_bond_line(o, n + 1, n + 2, d_sb, k_K_bond);
    n += 2;
  }
}

void output_ff_angle_line(std::ostream& o, int i, int j, int k, double a, double h)
{
  o << std::setw(8) << i << " " << std::setw(8) << j << " " << std::setw(8) << k << " "
    << std::setiosflags(std::ios_base::fixed) << std::setprecision(6) << std::setw(12)
    << a << " " << std::setprecision(1) << std::setw(8) << h << "\n";
}
void Chain::output_ffparm_angle(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double a = 0;
  Vec3d v1, v2;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-2; ++i) {
      v1 = v_residues_[i].get_cg_C_alpha().get_coordinate() - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      v2 = v_residues_[i+2].get_cg_C_alpha().get_coordinate() - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      a = vec_angle_deg (v1, v2);
      output_ff_angle_line(o, i + n + 1, i + n + 2, i + n + 3, a, k_K_angle);
    }
    n += n_residue_;
  } else {
    // ---------- angle BSP ----------
    v1 = v_residues_[0].get_cg_B().get_coordinate() - v_residues_[0].get_cg_S().get_coordinate();
    v2 = v_residues_[1].get_cg_P().get_coordinate() - v_residues_[0].get_cg_S().get_coordinate();
    a = vec_angle_deg (v1, v2);
    output_ff_angle_line(o, n + 2, n + 1, n + 3, a, k_K_angle);
    // ---------- angle SPS ----------
    v1 = v_residues_[0].get_cg_S().get_coordinate() - v_residues_[1].get_cg_P().get_coordinate();
    v2 = v_residues_[1].get_cg_S().get_coordinate() - v_residues_[1].get_cg_P().get_coordinate();
    a = vec_angle_deg (v1, v2);
    output_ff_angle_line(o, n + 1, n + 3, n + 4, a, k_K_angle);
    // -------------------- loop --------------------
    for (i = 1; i < n_residue_-1; ++i) {
      n += 3;
      // ---------- angle PSB ----------
      v1 = v_residues_[i].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i].get_cg_B().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      a = vec_angle_deg (v1, v2);
      output_ff_angle_line(o, n, n + 1, n + 2, a, k_K_angle);
      // ---------- angle PSP ----------
      v1 = v_residues_[i].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i+1].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      a = vec_angle_deg (v1, v2);
      output_ff_angle_line(o, n, n + 1, n + 3, a, k_K_angle);
      // ---------- angle BSP ----------
      v1 = v_residues_[i].get_cg_B().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i+1].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      a = vec_angle_deg (v1, v2);
      output_ff_angle_line(o, n + 2, n + 1, n + 3, a, k_K_angle);
      // ---------- angle SPS ----------
      v1 = v_residues_[i].get_cg_S().get_coordinate() - v_residues_[i+1].get_cg_P().get_coordinate();
      v2 = v_residues_[i+1].get_cg_S().get_coordinate() - v_residues_[i+1].get_cg_P().get_coordinate();
      a = vec_angle_deg (v1, v2);
      output_ff_angle_line(o, n + 1, n + 3, n + 4, a, k_K_angle);
    }
    n += 3;
    i = n_residue_ - 1;
    // ---------- angle PSB ----------
    v1 = v_residues_[i].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
    v2 = v_residues_[i].get_cg_B().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
    a = vec_angle_deg (v1, v2);
    output_ff_angle_line(o, n, n + 1, n + 2, a, k_K_angle);
    n += 2;
  }
}

void output_ff_dihedral_line(std::ostream& o, int i, int j, int k, int l, double d, double f, double g)
{
  o << std::setw(8) << i << " " << std::setw(8) << j << " "
    << std::setw(8) << k << " " << std::setw(8) << l << " "
    << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
    << std::setw(12) << d << " " << std::setprecision(1)
    << std::setw(8) << f << " " << std::setw(8) << g << "\n";
}
void Chain::output_ffparm_dihedral(std::ostream& o, int& n)
{
  if (chain_type_ == water || chain_type_ == other || chain_type_ == none)
    return;

  int i = 0;
  double d = 0;           // dihedral
  Vec3d v1, v2, v3, n1, n2;
  if (chain_type_ != DNA && chain_type_ != RNA && chain_type_ != na)
  {
    for (i = 0; i < n_residue_-3; ++i) {
      v1 = v_residues_[i].get_cg_C_alpha().get_coordinate() - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      v2 = v_residues_[i+2].get_cg_C_alpha().get_coordinate() - v_residues_[i+1].get_cg_C_alpha().get_coordinate();
      v3 = v_residues_[i+2].get_cg_C_alpha().get_coordinate() - v_residues_[i+3].get_cg_C_alpha().get_coordinate();
      n1 = v1 % v2;
      n2 = v2 % v3;
      d = vec_angle_deg (n1, n2);
      int ni = i + n;
      output_ff_dihedral_line(o, ni + 1, ni + 2, ni + 3, ni + 4, d, k_K_dihedral_1, k_K_dihedral_3);
    }
    n += n_residue_;
  } else {
    if (n_residue_ <= 2) {
      n += 5;
      return;
    }
    v1 = v_residues_[0].get_cg_S().get_coordinate() - v_residues_[1].get_cg_P().get_coordinate();
    v2 = v_residues_[1].get_cg_S().get_coordinate() - v_residues_[1].get_cg_P().get_coordinate();
    v3 = v_residues_[1].get_cg_S().get_coordinate() - v_residues_[2].get_cg_P().get_coordinate();
    n1 = v1 % v2;
    n2 = v2 % v3;
    d = vec_angle_deg (n1, n2);
    output_ff_dihedral_line(o, n + 1, n + 3, n + 4, n + 6, d, k_K_dihedral_1, k_K_dihedral_3);

    for (i = 1; i < n_residue_-1; ++i) {
      n += 3;
      // ---------- PSPS ----------
      v1 = v_residues_[i].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      v2 = v_residues_[i+1].get_cg_P().get_coordinate() - v_residues_[i].get_cg_S().get_coordinate();
      v3 = v_residues_[i+1].get_cg_P().get_coordinate() - v_residues_[i+1].get_cg_S().get_coordinate();
      n1 = v1 % v2;
      n2 = v2 % v3;
      d = vec_angle_deg (n1, n2);
      output_ff_dihedral_line(o, n, n + 1, n + 3, n + 4, d, k_K_dihedral_1, k_K_dihedral_3);

      if (i == n_residue_ - 2)
        break;
      // ---------- SPSP ----------
      v1 = v_residues_[i].get_cg_S().get_coordinate() - v_residues_[i+1].get_cg_P().get_coordinate();
      v2 = v_residues_[i+1].get_cg_S().get_coordinate() - v_residues_[i+1].get_cg_P().get_coordinate();
      v3 = v_residues_[i+1].get_cg_S().get_coordinate() - v_residues_[i+2].get_cg_P().get_coordinate();
      n1 = v1 % v2;
      n2 = v2 % v3;
      d = vec_angle_deg (n1, n2);
      output_ff_dihedral_line(o, n + 1, n + 3, n + 4, n + 6, d, k_K_dihedral_1, k_K_dihedral_3);
    }
    n += 5;
  }
}

void Chain::output_ffparm_native(std::ostream& o)
{
  int i = 0, j = 0;
  double d = -1, f = -1;
  ChainType cti, ctj;
  for (i = 0; i < n_residue_-4; ++i) {
    cti = v_residues_[i].get_chain_type();
    if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
      continue;
    for (j = i + 4; j < n_residue_; ++j) {
      ctj = v_residues_[j].get_chain_type();
      if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
        continue;
      d = residue_min_distance(v_residues_[i], v_residues_[j]);
      if ( d < g_cutoff)
      {
        f = residue_ca_distance(v_residues_[i], v_residues_[j]);
        o << std::setw(8) << i+1 << " " << std::setw(8) << j+1 << " "
          << std::setiosflags(std::ios_base::fixed) << std::setprecision(6)
          << std::setw(16) << f << " " << std::setprecision(5)
          << std::setw(12) << k_K_native << " " << "\n";
      }
    }
  }
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
    for (i = 0; i < s1; ++i) {
      c0.v_residues_.push_back(c1.v_residues_[i]);
      c0.n_residue_++;
    }
  if (s2 > 0)
    for (i = 0; i < s2; ++i) {
      c0.v_residues_.push_back(c2.v_residues_[i]);
      c0.n_residue_++;
    }
  return c0;
}

std::ostream& operator<<(std::ostream& o, Chain& c)
{
  int i = 0;
  int s = c.n_residue_;
  for (i = 0; i < s; ++i) {
    o << c.v_residues_[i];
  }
  o << "TER   " << "\n";
  return o;
}

}  // pinang

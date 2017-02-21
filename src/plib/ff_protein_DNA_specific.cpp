/*!
  @file ff_protein_DNA_specific.cpp
  @brief Basic class functions for PairProteinDNASpecific and PairProteinDNASpecificCombination.

  Basic class functions for PairProteinDNASpecific and
  PairProteinDNASpecificCombination, including read in amd simple calculations.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-08-01 11:49
  @copyright GNU Public License V3.0
*/

#include <fstream>
#include <sstream>
#include "ff_protein_DNA_specific.hpp"

namespace pinang {

PairProteinDNASpecific::PairProteinDNASpecific()
{
  protein_serial_     = -997;
  r_0_                = 9973.0;
  angle_0_0_          = -997.0;
  angle_53_0_         = -997.0;
  angle_NC_0_         = -997.0;
  sigma_              = 1.0;
  twice_sigma_square_ = 2.0;
  phi_                = 180.0;
  ene_pwm_A_          = 9973.0;
  ene_pwm_C_          = 9973.0;
  ene_pwm_G_          = 9973.0;
  ene_pwm_T_          = 9973.0;
}

std::istream& operator>>(std::istream& i, PairProteinDNASpecific& pairSS)
{
  double r0, angNC, ang0, ang53, sig, p, eneA, eneC, eneG, eneT;
  int serp;

  i >> serp >> r0 >> angNC >> ang0 >> ang53
    >> sig >> p >> eneA >> eneC >> eneG >> eneT;
  if (!i) return i;

  pairSS.protein_serial_     = serp - 1;
  pairSS.r_0_                = r0;
  pairSS.angle_NC_0_         = angNC;
  pairSS.angle_0_0_          = ang0;
  pairSS.angle_53_0_         = ang53;
  pairSS.sigma_              = sig;
  pairSS.twice_sigma_square_ = 2.0 * sig * sig;
  pairSS.phi_                = p;
  pairSS.ene_pwm_A_          = eneA * k_u_kB_T;
  pairSS.ene_pwm_C_          = eneC * k_u_kB_T;
  pairSS.ene_pwm_G_          = eneG * k_u_kB_T;
  pairSS.ene_pwm_T_          = eneT * k_u_kB_T;

  return i;
}

PairProteinDNASpecificCombination::PairProteinDNASpecificCombination()
{
  n_inter_pair_ = 0;
  protein_serial_ = -9973;
  interaction_pairs_.clear();
  d_cutoff_ = -1.0;
}

void PairProteinDNASpecificCombination::reset()
{
  n_inter_pair_ = 0;
  protein_serial_ = -9973;
  interaction_pairs_.clear();
  d_cutoff_ = -1.0;
}

int PairProteinDNASpecificCombination::add_interaction_pair(const PairProteinDNASpecific& pairSS)
{
  double tmp_cutoff = 0;
  if (n_inter_pair_ == 0) {
    interaction_pairs_.push_back(pairSS);
    protein_serial_ = pairSS.protein_serial_;
    n_inter_pair_++;
    tmp_cutoff = pairSS.r_0_ + 5 * pairSS.sigma_;
    d_cutoff_ = tmp_cutoff;
    return 0;
  } else {
    if (protein_serial_ == pairSS.protein_serial_) {
      interaction_pairs_.push_back(pairSS);
      tmp_cutoff = pairSS.r_0_ + 5 * pairSS.sigma_;
      d_cutoff_ = tmp_cutoff > d_cutoff_ ? tmp_cutoff : d_cutoff_;
      n_inter_pair_++;
      return 0;
    } else {
      return 1;
    }
  }
}

FFProteinDNASpecific::FFProteinDNASpecific()
{
  ss_pairwise_params_.clear();
  n_protein_particle_ = 0;
}

FFProteinDNASpecific::FFProteinDNASpecific(std::string ffp_file_name)
{
  ss_pairwise_params_.clear();
  std::ifstream ffp_file(ffp_file_name.c_str());
  std::string ffp_line;
  std::string tmp_s;
  PairProteinDNASpecific pair_tmp;
  PairProteinDNASpecificCombination pair_com_tmp;

  int tmp_i0;
  std::string::size_type m;

  while (ffp_file.good()) {
    std::getline(ffp_file, ffp_line);
    if (ffp_file.fail())
      break;

    std::istringstream tmp_sstr;
    m = ffp_line.find("[ protein-DNA seq-specific w/ PWM ]");
    if (m != std::string::npos) {
      tmp_sstr.str(ffp_line);
      for (int i = 0; i < 6; ++i)
        tmp_sstr >> tmp_s;
      tmp_sstr >> tmp_i0;  // read in total number of pro-DNA interactions;
      getline(ffp_file, ffp_line);  // read in a comment line;
      for (int i = 0; i < tmp_i0; ++i) {  // read in all interaction details...;
        ffp_file >> pair_tmp;
        if (pair_com_tmp.add_interaction_pair(pair_tmp)) {
          if (pair_com_tmp.n_inter_pair_ > 0) {
            ss_pairwise_params_.push_back(pair_com_tmp);
          }
          pair_com_tmp.reset();
          pair_com_tmp.add_interaction_pair(pair_tmp);
        }
      }
      if (pair_com_tmp.n_inter_pair_ > 0) {
        ss_pairwise_params_.push_back(pair_com_tmp);
      }
      break;
    } else {
      m = ffp_line.find("[ protein-DNA seq-specific ]");
      if (m != std::string::npos) {
        std::cout << " Please apply patch of PWM to protein-DNA seq-specific interactions ! "
                  << "\n";
        exit(EXIT_SUCCESS);
      }
    }
  }
  ffp_file.close();

  n_protein_particle_ = ss_pairwise_params_.size();
}

}  // pinang

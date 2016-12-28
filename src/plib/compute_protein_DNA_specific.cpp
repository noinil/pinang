/*!
  @file compute_protein_DNA_specific.cpp
  @brief Compute energy and force for protein-DNA sequence specific interactions.

  Compute protein-DNA sequence specific interaction energy and forces.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-08-01 15:54
  @copyright GNU Public License V3.0
*/

#include <cmath>
#include <iomanip>
#include "ff_protein_DNA_specific.hpp"

namespace pinang {

double FFProteinDNASpecific::compute_energy_protein_DNA_specific(Topology &top, Conformation &conf)
{
  double total_energy = 0;
  std::vector<int> dna_index;

  int i, j, k;
  if (conf.get_size() != top.get_size()) {
    std::cout << " Inconsistent particle number in Topology and Conformation" << "\n";
    exit(EXIT_SUCCESS);
  }
  int n_parts = top.get_size();
  for (i = 0; i < n_parts; ++i) {
    if (top.get_particle(i).get_atom_name() == "DB  ") {
      dna_index.push_back(i);
    }
  }

  // std::cout << "n_protein_particle_ :  " << n_protein_particle_ << "\n";
  for (i = 0; i < n_protein_particle_; ++i) {
    PairProteinDNASpecificCombination tmp_pair = ss_pairwise_params_[i];
    int proi = tmp_pair.protein_serial_;  // protein particle index;
    // std::cout << " protein particle: " << i << "   serial: " << proi + 1 << "\n";
    Vec3d tmp_c_CA = conf.get_coordinate(proi);         // protein Calpha coordinates;
    Vec3d tmp_c_CA_N, tmp_c_CA_C;
    int calpha_term_N, calpha_term_C;
    k = proi - 1;
    if (k < 0 || top.get_particle(k).get_chain_ID() != top.get_particle(proi).get_chain_ID()) {
      tmp_c_CA_N = tmp_c_CA;
      calpha_term_N = 1;
    } else {
      calpha_term_N = 0;
      tmp_c_CA_N = conf.get_coordinate(k);  // Coor of N' CA
    }
    k = proi + 1;
    if (k >= top.get_size() || top.get_particle(k).get_chain_ID() != top.get_particle(proi).get_chain_ID()) {
      tmp_c_CA_C = tmp_c_CA;
      calpha_term_C = 1;
    } else {
      calpha_term_C = 0;
      tmp_c_CA_C = conf.get_coordinate(k);  // Coor of C' CA
    }
    if (calpha_term_N * calpha_term_C > 0) {
      std::cout << " Single residue Chain!!! WTF!!! \n";
      exit(EXIT_SUCCESS);
    }
    Vec3d tmp_CCA_NCA = tmp_c_CA_N - tmp_c_CA_C;
    double dist_cutoff = tmp_pair.d_cutoff_;
    for (j = 0; j < dna_index.size(); ++j) {
      int dnai = dna_index[j];                          // DNA particle index;
      char tmp_chain_id = top.get_particle(dnai).get_chain_ID();
      std::string tmp_base_name = top.get_particle(dnai).get_residue_name();
      Vec3d tmp_c_B0 = conf.get_coordinate(dnai);       // DNA Base coordinates;
      Vec3d tmp_B0_CA = tmp_c_CA - tmp_c_B0;            // vector from Base to Calpha;
      double tmp_distance = tmp_B0_CA.norm();
      if (tmp_distance >= dist_cutoff) {
        continue;
      }
      Vec3d tmp_c_S0 = conf.get_coordinate(dnai - 1);
      double tmp_angle_0 = vec_angle_deg(tmp_c_S0 - tmp_c_B0, tmp_B0_CA);
      Vec3d tmp_c_B5, tmp_c_B3;
      if (dnai - 3 < 0 || top.get_particle(dnai - 3).get_chain_ID() != tmp_chain_id) {
        tmp_c_B5 = conf.get_coordinate(dnai);
      } else {
        tmp_c_B5 = conf.get_coordinate(dnai - 3);
      }
      if (dnai + 3 >= top.get_size() || top.get_particle(dnai + 3).get_chain_ID() != tmp_chain_id) {
        tmp_c_B3 = conf.get_coordinate(dnai);
      } else {
        tmp_c_B3 = conf.get_coordinate(dnai + 3);
      }
      Vec3d tmp_B5_B3 = tmp_c_B3 - tmp_c_B5;
      double tmp_dist_B5_B3 = tmp_B5_B3.norm();
      if (tmp_dist_B5_B3 < 0.0001) {
        continue;
      }
      double tmp_angle_53 = vec_angle_deg(tmp_B5_B3, tmp_B0_CA);
      double tmp_angle_NC = vec_angle_deg(tmp_CCA_NCA,  tmp_B0_CA);
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CORE CALCULATION! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      for (k = 0; k < tmp_pair.n_inter_pair_; ++k) {
        PairProteinDNASpecific p = tmp_pair.interaction_pairs_[k];
        double f1 = 0, f2 = 0, f3 = 0, f4 = 0;  // f1: bond; f2: angle 0; f3: angle NC; f4: angle 53;
        double dr = tmp_distance - p.r_0_;
        f1 = exp(- (dr * dr) / p.twice_sigma_square_);
        double delta_theta_0 = std::abs(tmp_angle_0 - p.angle_0_0_);
        double delta_theta_NC = std::abs(tmp_angle_NC - p.angle_NC_0_);
        double delta_theta_53 = std::abs(tmp_angle_53 - p.angle_53_0_);
        if (delta_theta_0 < p.phi_) {
          f2 = 1;
        } else if (delta_theta_0 < p.phi_ * 2.0) {
          double cos_theta_0 = cos(delta_theta_0 / 180.0 * 3.14159265);
          f2 = 1 - cos_theta_0 * cos_theta_0;
        } else {
          continue;
        }
        if (delta_theta_NC < p.phi_) {
          f3 = 1;
        } else if (delta_theta_NC < p.phi_ * 2.0) {
          double cos_theta_NC = cos(delta_theta_NC / 180.0 * 3.14159265);
          f3 = 1 - cos_theta_NC * cos_theta_NC;
        } else {
          continue;
        }
        if (delta_theta_53 < p.phi_) {
          f4 = 1;
        } else if (delta_theta_53 < p.phi_ * 2.0) {
          double cos_theta_53 = cos(delta_theta_53 / 180.0 * 3.14159265);
          f4 = 1 - cos_theta_53 * cos_theta_53;
        } else {
          continue;
        }
        double f = f1 * f2 * f3 * f4;
        double ene_pwm_base = 0;
        if (tmp_base_name == "DA ") {
          ene_pwm_base = p.ene_pwm_A_;
        } else if (tmp_base_name == "DC ") {
          ene_pwm_base = p.ene_pwm_C_;
        } else if (tmp_base_name == "DG ") {
          ene_pwm_base = p.ene_pwm_G_;
        } else if (tmp_base_name == "DT ") {
          ene_pwm_base = p.ene_pwm_T_;
        }
        double e = ene_pwm_base * f;
        total_energy += e;
      }
    }
  }
  // std::cout << " ============~~~~~~~~~~~~~~~==============" << "\n";
  // std::cout << total_energy << "\n";
  // std::cout << " ============~~~~~~~~~~~~~~~==============" << "\n";

  return total_energy;
}

}  // pinang

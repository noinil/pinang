/*!
  @file compute_protein_DNA_specific.cpp
  @brief Compute energy and force for protein-DNA sequence specific interactions.

  Compute protein-DNA sequence specific interaction energy and forces.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-08-01 15:54
  @copyright GNU Public License V3.0
*/

#include <cmath>
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

  for (i = 0; i < n_protein_particle_; ++i) {
    PairProteinDNASpecificCombination tmp_pair = ss_pairwise_params_[i];
    int proi = tmp_pair.protein_serial_;  // protein particle index;
    std::cout << " protein particle: " << i << "   serial: " << proi + 1 << "\n";
    Vec3d tmp_c_CA = conf.get_coordinate(proi);         // protein Calpha coordinates;
    double dist_cutoff = tmp_pair.d_cutoff_;
    for (j = 0; j < dna_index.size(); ++j) {
      int dnai = dna_index[j];                          // DNA particle index;
      char tmp_chain_id = top.get_particle(dnai).get_chain_ID();
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
      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CORE CALCULATION! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double tmp_angle_factor = 0;
      double tmp_factor = 0;
      double best_d = 0, best_a0 = 0, best_a53 = 0;
      for (k = 0; k < tmp_pair.n_inter_pair_; ++k) {
        PairProteinDNASpecific p = tmp_pair.interaction_pairs_[k];
        double f1, f2, f3;
        double dr = tmp_distance - p.r_0_;
        f1 = exp(- (dr * dr) / p.twice_sigma_square_);
        double delta_theta_0 = std::abs(tmp_angle_0 - p.angle_0_0_);
        double delta_theta_53 = std::abs(tmp_angle_53 - p.angle_53_0_);
        if (delta_theta_0 < p.phi_) {
          f2 = 1;
        } else if (delta_theta_0 < p.phi_ * 2.0) {
          double cos_theta_0 = cos(delta_theta_0 / 180.0 * 3.14159265);
          f2 = 1 - cos_theta_0 * cos_theta_0;
        } else {
          f2 = 0;
        }
        if (delta_theta_53 < p.phi_) {
          f3 = 1;
        } else if (delta_theta_53 < p.phi_ * 2.0) {
          double cos_theta_53 = cos(delta_theta_53 / 180.0 * 3.14159265);
          f3 = 1 - cos_theta_53 * cos_theta_53;
        } else {
          f3 = 0;
        }
        double f = f1 * f2 * f3;
        double fangle = f2 * f3;
        if (fangle > tmp_angle_factor) {
          tmp_factor = f;
          tmp_angle_factor = fangle;
          best_d = tmp_distance;
          best_a0 = tmp_angle_0;
          best_a53 = tmp_angle_53;
        }
      }
      double tmp_energy = -1.0 * tmp_factor;
      total_energy += tmp_energy;
      std::cout << "   " << proi + 1 << " - " << dnai + 1 << "   " << tmp_energy
                << "     d: " << best_d << "      a0: " << best_a0 << "      a53: " << best_a53
                << "\n";
    }
  }
  std::cout << " ============~~~~~~~~~~~~~~~==============" << "\n";
  std::cout << total_energy << "\n";

  return total_energy;
}

}  // pinang

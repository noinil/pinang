/*!
  @file ff_protein_DNA_specific.hpp
  @brief Protein-DNA sequence specific interactions.

  In this file we define a class that describe force field parameter and formula
  of protein-DNA sequence specific interactions.
  
  @author Cheng Tan (noinil@gmail.com)
  @date 2016-08-01 11:07
  @copyright GNU Public License V3.0
*/

#include "geometry.hpp"
#include "topology.hpp"

#ifndef PINANG_FF_PROTEIN_DNA_SPECIFIC_H
#define PINANG_FF_PROTEIN_DNA_SPECIFIC_H

namespace pinang {

class PairProteinDNASpecific;
class PairProteinDNASpecificCombination;

/*!
  @brief Force field details of protein-DNA sequence-specific interactions.

  This class defines formula and parameters for protein-DNA sequence-specific
  interactions.
*/
class FFProteinDNASpecific
{
 public:
  //! @brief Create empty FFProteinDNASpecific object.
  //! @retval A FFProteinDNASpecific object.
  FFProteinDNASpecific();
  //! @brief Create FFProteinDNASpecific object from .ffp file.
  //! @retval A FFProteinDNASpecific object with parameters set.
  FFProteinDNASpecific(std::string);
  virtual ~FFProteinDNASpecific() {ss_pairwise_params_.clear();}

  //! @brief Compute protein-DNA sequence specific interaction energy.
  //! @param Topology and conformation.
  //! @retval Energy.
  double compute_energy_protein_DNA_specific(Topology&, Conformation&);
 protected:
  int n_protein_particle_;
  std::vector<PairProteinDNASpecificCombination> ss_pairwise_params_;
};

/*!
  @brief PairProteinDNASpecificCombination

  Store several protein-DNA interaction pairs for one protein particle.
*/
class PairProteinDNASpecificCombination
{
 public:
  //! @brief Create empty PairProteinDNASpecific object.
  //! @retval A PairProteinDNASpecific object.
  PairProteinDNASpecificCombination();
  virtual ~PairProteinDNASpecificCombination() {interaction_pairs_.clear();}

  //! @brief Add a new interaction pair.
  //! @param PairProteinDNASpecific.
  //! @return Status of adding PairProteinDNASpecific.
  //! @retval 0: Success.
  int add_interaction_pair(const PairProteinDNASpecific&);

  //! @brief Reset properties of PairProteinDNASpecificCombination.
  void reset();

  friend class FFProteinDNASpecific;
 protected:
  std::vector<PairProteinDNASpecific> interaction_pairs_;  //!< pro-DNA interaction pairs for one pro particle.
  int n_inter_pair_;                                      //!< number of pro-DNA inter pairs.
  int protein_serial_;  //!< protein particle serial number.
  double d_cutoff_;     //!< protein-DNA distance cutoff (= r0 + 5.0 )
};

/*!
  @brief "Pairwise" parameters for protein-DNA sequence specific interactions.

  Store protein index, r_0, angle_0, angle_53, sigma, and phi.
*/
class PairProteinDNASpecific
{
 public:
  //! @brief Create empty PairProteinDNASpecific object.
  //! @retval A PairProteinDNASpecific object.
  PairProteinDNASpecific();
  virtual ~PairProteinDNASpecific() {};

  //! @brief Read in PairProteinDNASpecific parameters.
  friend std::istream& operator>>(std::istream&, PairProteinDNASpecific&);
  friend class PairProteinDNASpecificCombination;
  friend class FFProteinDNASpecific;

 protected:
  int    protein_serial_;               //!< protein particle serial number.
  double r_0_;                          //!< distance in the native structure.
  double angle_0_0_;                    //!< angle sugar-base-Calpha in the native structure.
  double angle_53_0_;                   //!< angle (5'base-3'base)--(base-Calpha) in the native structure.
  double angle_NC_0_;                   //!< angle (N'CA-C'CA)--(Calpha-base) in the native structure.
  double sigma_;                        //!< sigma in the Gaussian potential function of distance.
  double twice_sigma_square_;           //!< 2 * sigma * sigma.
  double phi_;                          //!< angle threshold in the potential of angles.

  double ene_pwm_A_;                    //!< PWM energy for base A.
  double ene_pwm_C_;                    //!< PWM energy for base C.
  double ene_pwm_G_;                    //!< PWM energy for base G.
  double ene_pwm_T_;                    //!< PWM energy for base T.
};

}  // pinang

#endif

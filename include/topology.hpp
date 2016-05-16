/*!
************************************************************
@file topology.hpp
@brief Definition of class Topology.

In this file class Topology is defined.  Topology is nothing more than a
collection of particles.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:11
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_TOPOLOGY_H
#define PINANG_TOPOLOGY_H

#include "particle.hpp"

#include <vector>
#include <string>

namespace pinang {

/*!
************************************************************
@brief A set of particles.

The class Topology only contains physical properties of a set of particles.

@todo Add interaction (such as bonds, angles, dihedrals...) information.
************************************************************
*/

class Topology
{
 public:
  Topology(const std::string& s);
  virtual ~Topology() {particles_.clear();}

  void reset();

  int get_size() { return n_particle_; }
  Particle& get_particle(int);

 protected:
  std::vector<Particle> particles_;  //!< A set of Particle objects.
  int n_particle_;                   //!< Number of particles in Topology.
};

}

#endif

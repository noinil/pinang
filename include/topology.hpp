// -*-c++-*-

#ifndef PINANG_TOPOLOGY_H
#define PINANG_TOPOLOGY_H

#include "constants.hpp"
#include "particle.hpp"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

namespace pinang {

class Topology
{
 public:
  Topology(const std::string& s);
  virtual ~Topology() {particles_.clear();}

  void reset();

  int get_size() { return n_particle_; }
  Particle& get_particle(int);

 protected:
  std::vector<Particle> particles_;
  int n_particle_;
};

}

#endif

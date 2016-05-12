// -*-c++-*-

#ifndef PINANG_TOPOLOGY_H
#define PINANG_TOPOLOGY_H

#include "constants.h"
#include "particle.h"

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

  inline void reset();

  inline int get_size() {return n_particle_;}
  Particle& get_particle(int);

 protected:
  std::vector<Particle> particles_;
  int n_particle_;
};

Topology::Topology(const std::string& s)
{
  n_particle_ = 0;
  particles_.clear();

  Particle p;

  std::ifstream ifile(s.c_str());
  if (!ifile.is_open())
  {
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cerr << " ERROR: Cannot read top file: " << s << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string inp_line;
  while (ifile.good()) {
    std::getline(ifile, inp_line);

    if (ifile.fail())
      break;

    std::size_t found = inp_line.find("particles");
    if (found!=std::string::npos){
      std::string stmp;
      std::istringstream tmp_sstr;
      tmp_sstr.str ( inp_line );
      tmp_sstr >> stmp  >> stmp  >> stmp
               >> n_particle_;
      std::getline(ifile, inp_line);
      for (int i = 0; i < n_particle_ ; i++) {
        ifile >> p;
        particles_.push_back(p);
      }
      std::cout << " Total particle number: "
                << n_particle_
                << " in top file: " << s
                << std::endl;
      break;
    }
  }
  ifile.close();
}

inline void Topology::reset()
{
  n_particle_ = 0;
  particles_.clear();
}

Particle& Topology::get_particle(int n)
{
  if ( n >= n_particle_ || n < 0)
  {
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cerr << " ERROR: Atom index out of range in Topology. " << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    return particles_[n];
  }
}

}
#endif

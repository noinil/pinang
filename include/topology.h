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
  virtual ~Topology() {_particles.clear();}

  inline void reset();

  int get_size() {return _n_particle;}
  const Particle& particle(int n) const;

 protected:
  std::vector<Particle> _particles;
  int _n_particle;
};

Topology::Topology(const std::string& s)
{
  _n_particle = 0;
  _particles.clear();

  Particle p;

  std::ifstream ifile(s.c_str());
  if (!ifile.is_open())
  {
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
               >> _n_particle;
      std::getline(ifile, inp_line);
      for (int i = 0; i < _n_particle ; i++) {
        ifile >> p;
        _particles.push_back(p);
      }
      std::cout << " Total particle number: "
                << _n_particle
                << " in top file: " << s
                << std::endl;
      break;
    }
  }
  ifile.close();
}

inline void Topology::reset()
{
  _n_particle = 0;
  _particles.clear();
}

const Particle& Topology::particle(int n) const
{
  if ( n >= _n_particle || n < 0)
  {
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << std::endl;
    std::cout << " ============================================== " << std::endl;
    std::cerr << " ERROR: Atom index out of range in Topology. " << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    return _particles[n];
  }
}

}
#endif

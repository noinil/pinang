/*!
@file topology.cpp
@brief Define functions of class Topology.

Definitions of member or friend functions of class Topology.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-24 15:46
@copyright GNU Public License V3.0
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include "topology.hpp"

namespace pinang {

Topology::Topology(const std::string& s)
{
  n_particle_ = 0;
  v_particles_.clear();

  Particle p;

  std::ifstream ifile(s.c_str());
  if (!ifile.is_open())
  {
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << "\n";
    std::cerr << " ERROR: Cannot read top file: " << s << "\n";
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
        v_particles_.push_back(p);
      }
      std::cout << " Total particle number: "
                << n_particle_
                << " in top file: " << s
                << "\n";
      break;
    }
  }
  ifile.close();
}

void Topology::reset()
{
  n_particle_ = 0;
  v_particles_.clear();
}

Particle& Topology::get_particle(int n)
{
  if ( n >= n_particle_ || n < 0)
  {
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << "\n";
    std::cerr << " ERROR: Atom index out of range in Topology. " << "\n";
    exit(EXIT_SUCCESS);
  } else {
    return v_particles_[n];
  }
}

}  // pinang

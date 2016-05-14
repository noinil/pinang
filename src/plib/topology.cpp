#include <iostream>
#include <sstream>
#include <fstream>
#include "topology.hpp"

namespace pinang {

Topology::Topology(const std::string& s)
{
  n_particle_ = 0;
  particles_.clear();

  Particle p;

  std::ifstream ifile(s.c_str());
  if (!ifile.is_open())
  {
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << std::endl;
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

void Topology::reset()
{
  n_particle_ = 0;
  particles_.clear();
}

Particle& Topology::get_particle(int n)
{
  if ( n >= n_particle_ || n < 0)
  {
    std::cout << " ~             PINANG :: TOPOLOGY             ~ " << std::endl;
    std::cerr << " ERROR: Atom index out of range in Topology. " << std::endl;
    exit(EXIT_SUCCESS);
  } else {
    return particles_[n];
  }
}

}  // pinang

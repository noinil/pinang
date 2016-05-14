// -*-c++-*-

#ifndef PINANG_CONSTANTS_H_
#define PINANG_CONSTANTS_H_

#include <map>
#include <string>
#include <iostream>
// #include <cstdlib>

namespace pinang {

extern double g_cutoff;

const long double k_pi = 3.14159265358979323846;

// Units.
const double k_u_mass = 1.0;

// CG MD Parameters.
const double k_K_bond = 100.0;
const double k_K_angle = 20.0;
const double k_K_dihedral_1 = 1.0;
const double k_K_dihedral_3 = 0.5;
const double k_K_native = 1.0;
const double k_K_nonnative = 1.0;

enum ChainType {none=0, protein=1, DNA=2, RNA=3, water=4, ion=5, other=6, na=7};

class PhysicalProperty
{
 public:
  PhysicalProperty();
  ~PhysicalProperty();

  std::string get_short_name(const std::string&);
  double get_charge(const std::string&);
  double get_mass(const std::string&);
  ChainType get_chain_type(const std::string&);

 private:
  std::map<std::string, std::string> map_resName_shortName;
  std::map<std::string, double> map_resName_charge;
  std::map<std::string, double> map_resName_mass;
  std::map<std::string, ChainType> map_resName_chainType;
};

}

#endif

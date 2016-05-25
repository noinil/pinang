/*!
************************************************************
@file constants.cpp
@brief Defines lots of physical constants.

Definitions of physical constants and functions.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-24 15:39
@copyright GNU Public License V3.0
************************************************************
*/

#include "constants.hpp"

namespace pinang {

double g_cutoff = 6.5;

PhysicalProperty::PhysicalProperty()
{
  // ---------- full name to short name ----------
  map_resName_shortName["ALA"] = "aA";
  map_resName_shortName["ARG"] = "aR";
  map_resName_shortName["ASN"] = "aN";
  map_resName_shortName["ASP"] = "aD";
  map_resName_shortName["CYS"] = "aC";
  map_resName_shortName["GLN"] = "aQ";
  map_resName_shortName["GLU"] = "aE";
  map_resName_shortName["GLY"] = "aG";
  map_resName_shortName["HIS"] = "aH";
  map_resName_shortName["ILE"] = "aI";
  map_resName_shortName["LEU"] = "aL";
  map_resName_shortName["LYS"] = "aK";
  map_resName_shortName["MET"] = "aM";
  map_resName_shortName["PHE"] = "aF";
  map_resName_shortName["PRO"] = "aP";
  map_resName_shortName["SER"] = "aS";
  map_resName_shortName["SEC"] = "aU";
  map_resName_shortName["THR"] = "aT";
  map_resName_shortName["TRP"] = "aW";
  map_resName_shortName["TYR"] = "aY";
  map_resName_shortName["VAL"] = "aV";
  map_resName_shortName["A  "] = "nA";
  map_resName_shortName["DA "] = "dA";
  map_resName_shortName["DA3"] = "dA";
  map_resName_shortName["DA5"] = "dA";
  map_resName_shortName["G  "] = "nG";
  map_resName_shortName["DG "] = "dG";
  map_resName_shortName["DG3"] = "dG";
  map_resName_shortName["DG5"] = "dG";
  map_resName_shortName["C  "] = "nC";
  map_resName_shortName["DC "] = "dC";
  map_resName_shortName["DC3"] = "dC";
  map_resName_shortName["DC5"] = "dC";
  map_resName_shortName["T  "] = "nT";
  map_resName_shortName["DT "] = "dT";
  map_resName_shortName["DT3"] = "dT";
  map_resName_shortName["DT5"] = "dT";
  map_resName_shortName["RA "] = "rA";
  map_resName_shortName["RG "] = "rG";
  map_resName_shortName["RC "] = "rC";
  map_resName_shortName["RU "] = "rU";
  map_resName_shortName["U  "] = "nU";
  map_resName_shortName["ZN "] = "iz";
  map_resName_shortName["CA "] = "ic";
  map_resName_shortName["MG "] = "im";
  map_resName_shortName["HOH"] = "wt";

  // ---------- full name to mass ----------
  map_resName_mass["ALA"] = 71.0788;
  map_resName_mass["ARG"] = 156.1875;
  map_resName_mass["ASN"] = 114.1038;
  map_resName_mass["ASP"] = 115.0886;
  map_resName_mass["CYS"] = 103.1388;
  map_resName_mass["GLU"] = 129.1155;
  map_resName_mass["GLN"] = 128.1307;
  map_resName_mass["GLY"] = 57.0519;
  map_resName_mass["HIS"] = 137.1411;
  map_resName_mass["ILE"] = 113.1594;
  map_resName_mass["LEU"] = 113.1594;
  map_resName_mass["LYS"] = 128.1741;
  map_resName_mass["MET"] = 131.1926;
  map_resName_mass["PHE"] = 147.1766;
  map_resName_mass["PRO"] = 97.1167;
  map_resName_mass["SER"] = 87.0782;
  map_resName_mass["SEC"] = 150.0379;
  map_resName_mass["THR"] = 101.1051;
  map_resName_mass["TRP"] = 186.2132;
  map_resName_mass["TYR"] = 163.1760;
  map_resName_mass["VAL"] = 99.1326;
  map_resName_mass["A  "] = 312.2;
  map_resName_mass["DA "] = 312.2;
  map_resName_mass["DA3"] = 312.2;
  map_resName_mass["DA5"] = 312.2;
  map_resName_mass["C  "] = 288.17;
  map_resName_mass["DC "] = 288.17;
  map_resName_mass["DC3"] = 288.17;
  map_resName_mass["DC5"] = 288.17;
  map_resName_mass["G  "] = 328.2;
  map_resName_mass["DG "] = 328.2;
  map_resName_mass["DG3"] = 328.2;
  map_resName_mass["DG5"] = 328.2;
  map_resName_mass["T  "] = 303.171;
  map_resName_mass["DT "] = 303.171;
  map_resName_mass["DT3"] = 303.171;
  map_resName_mass["DT5"] = 303.171;
  map_resName_mass["RA "] = 328.198;
  map_resName_mass["RC "] = 304.173;
  map_resName_mass["RU "] = 305.158;
  map_resName_mass["RG "] = 344.197;
  map_resName_mass["U  "] = 305.158;
  map_resName_mass["ZN "] = 65.409;
  map_resName_mass["CA "] = 40.08;
  map_resName_mass["MG "] = 24.305;
  map_resName_mass["HOH"] = 18.014;

  // ---------- full name to charg ----------
  map_resName_charge["ALA"] = 0.0;
  map_resName_charge["ARG"] = 1.0;
  map_resName_charge["ASN"] = 0.0;
  map_resName_charge["ASP"] = -1.0;
  map_resName_charge["CYS"] = 0.0;
  map_resName_charge["GLU"] = -1.0;
  map_resName_charge["GLN"] = 0.0;
  map_resName_charge["GLY"] = 0.0;
  map_resName_charge["HIS"] = 1.0;
  map_resName_charge["ILE"] = 0.0;
  map_resName_charge["LEU"] = 0.0;
  map_resName_charge["LYS"] = 1.0;
  map_resName_charge["MET"] = 0.0;
  map_resName_charge["PHE"] = 0.0;
  map_resName_charge["PRO"] = 0.0;
  map_resName_charge["SER"] = 0.0;
  map_resName_charge["SEC"] = 0.0;
  map_resName_charge["THR"] = 0.0;
  map_resName_charge["TRP"] = 0.0;
  map_resName_charge["TYR"] = 0.0;
  map_resName_charge["VAL"] = 0.0;
  map_resName_charge["A  "] = -1.0;
  map_resName_charge["DA "] = -1.0;
  map_resName_charge["DA3"] = -1.0;
  map_resName_charge["DA5"] = -1.0;
  map_resName_charge["C  "] = -1.0;
  map_resName_charge["DC "] = -1.0;
  map_resName_charge["DC3"] = -1.0;
  map_resName_charge["DC5"] = -1.0;
  map_resName_charge["G  "] = -1.0;
  map_resName_charge["DG "] = -1.0;
  map_resName_charge["DG3"] = -1.0;
  map_resName_charge["DG5"] = -1.0;
  map_resName_charge["T  "] = -1.0;
  map_resName_charge["DT "] = -1.0;
  map_resName_charge["DT3"] = -1.0;
  map_resName_charge["DT5"] = -1.0;
  map_resName_charge["RA "] = -1.0;
  map_resName_charge["RC "] = -1.0;
  map_resName_charge["RG "] = -1.0;
  map_resName_charge["RU "] = -1.0;
  map_resName_charge["U  "] = -1.0;
  map_resName_charge["ZN "] = 2.0;
  map_resName_charge["CA "] = 2.0;
  map_resName_charge["MG "] = 2.0;
  map_resName_charge["HOH"] = 0.0;

  // ---------- full name to residue type ----------
  map_resName_chainType["ALA"] = protein;
  map_resName_chainType["ARG"] = protein;
  map_resName_chainType["ASN"] = protein;
  map_resName_chainType["ASP"] = protein;
  map_resName_chainType["CYS"] = protein;
  map_resName_chainType["GLU"] = protein;
  map_resName_chainType["GLN"] = protein;
  map_resName_chainType["GLY"] = protein;
  map_resName_chainType["HIS"] = protein;
  map_resName_chainType["ILE"] = protein;
  map_resName_chainType["LEU"] = protein;
  map_resName_chainType["LYS"] = protein;
  map_resName_chainType["MET"] = protein;
  map_resName_chainType["PHE"] = protein;
  map_resName_chainType["PRO"] = protein;
  map_resName_chainType["SER"] = protein;
  map_resName_chainType["SEC"] = protein;
  map_resName_chainType["THR"] = protein;
  map_resName_chainType["TRP"] = protein;
  map_resName_chainType["TYR"] = protein;
  map_resName_chainType["VAL"] = protein;
  map_resName_chainType["A  "] = na;
  map_resName_chainType["DA "] = DNA;
  map_resName_chainType["DA3"] = DNA;
  map_resName_chainType["DA5"] = DNA;
  map_resName_chainType["T  "] = na;
  map_resName_chainType["DT "] = DNA;
  map_resName_chainType["DT3"] = DNA;
  map_resName_chainType["DT5"] = DNA;
  map_resName_chainType["G  "] = na;
  map_resName_chainType["DG "] = DNA;
  map_resName_chainType["DG3"] = DNA;
  map_resName_chainType["DG5"] = DNA;
  map_resName_chainType["C  "] = na;
  map_resName_chainType["DC "] = DNA;
  map_resName_chainType["DC3"] = DNA;
  map_resName_chainType["DC5"] = DNA;
  map_resName_chainType["RA "] = RNA;
  map_resName_chainType["RC "] = RNA;
  map_resName_chainType["RG "] = RNA;
  map_resName_chainType["RU "] = RNA;
  map_resName_chainType["U  "] = na;
  map_resName_chainType["ZN "] = ion;
  map_resName_chainType["CA "] = ion;
  map_resName_chainType["MG "] = ion;
  map_resName_chainType["HOH"] = water;
}

PhysicalProperty::~PhysicalProperty() {
  map_resName_shortName.clear();
  map_resName_charge.clear();
  map_resName_mass.clear();
  map_resName_chainType.clear();
}

std::string PhysicalProperty::get_short_name(const std::string& s)
{
  auto search = map_resName_shortName.find(s);
  if (search != map_resName_shortName.end()) {
    return search->second;
  } else {
    std::cout << " ~             PINANG :: constants.hpp          ~ " << "\n";
    std::cerr << " ERROR: Cannot find residue name for short name: "
              << s << "\n";
    exit(EXIT_FAILURE);
  }
}

double PhysicalProperty::get_charge(const std::string& s)
{
  auto search = map_resName_charge.find(s);
  if (search != map_resName_charge.end()) {
    return search->second;
  } else {
    std::cout << " ~             PINANG :: constants.hpp          ~ " << "\n";
    std::cerr << " ERROR: Cannot find residue name for charge: "
              << s << "\n";
    exit(EXIT_FAILURE);
  }
}

double PhysicalProperty::get_mass(const std::string& s)
{
  auto search = map_resName_mass.find(s);
  if (search != map_resName_mass.end()) {
    return search->second;
  } else {
    std::cout << " ~             PINANG :: constants.hpp          ~ " << "\n";
    std::cerr << " ERROR: Cannot find residue name for mass: "
              << s << "\n";
    exit(EXIT_FAILURE);
  }
}

ChainType PhysicalProperty::get_chain_type(const std::string& s)
{
  auto search = map_resName_chainType.find(s);
  if (search != map_resName_chainType.end()) {
    return search->second;
  } else {
    std::cout << " ~             PINANG :: constants.hpp          ~ " << "\n";
    std::cerr << " ERROR: Cannot find residue name for chain type: "
              << s << "\n";
    exit(EXIT_FAILURE);
  }
}


}  // pinang

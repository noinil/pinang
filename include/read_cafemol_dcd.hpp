#ifndef PINANG_READ_DCD_H_
#define PINANG_READ_DCD_H_

#include "conformation.hpp"

#include <vector>
#include <fstream>

namespace pinang {

int read_cafemol_dcd(std::ifstream&, std::vector<Conformation>&);

}
#endif

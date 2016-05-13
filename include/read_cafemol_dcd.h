// -*-c++-*-

#ifndef PINANG_READ_DCD_H_
#define PINANG_READ_DCD_H_

#include "conformation.h"

#include <vector>
#include <fstream>

namespace pinang {

int read_cafemol_dcd(std::ifstream&, std::vector<Conformation>&);

}
#endif

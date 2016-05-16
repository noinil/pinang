/*!
************************************************************
@file read_cafemol_dcd.hpp
@brief Basic function of reading CafeMol style dcd files.

In this file a function that can read CafeMol dcd file is provided.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-16 18:08
@copyright GNU Public License V3.0
************************************************************
*/

#ifndef PINANG_READ_DCD_H_
#define PINANG_READ_DCD_H_

#include "conformation.hpp"

#include <vector>
#include <fstream>

namespace pinang {

int read_cafemol_dcd(std::ifstream&, std::vector<Conformation>&);

}
#endif

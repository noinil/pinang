/*!
  @file read_cafemol_dcd.hpp
  @brief Basic function of reading CafeMol style dcd files.

  In this file a function that can read CafeMol dcd file is provided.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-16 18:08
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_READ_DCD_H_
#define PINANG_READ_DCD_H_

#include "conformation.hpp"

#include <vector>
#include <fstream>

namespace pinang {

//! @brief Read DCD information into Conformations.
//! @param DCD file.
//! @param Conformation.
//! @return Status of reading DCD file.
//! @retval 1: Failure.
//! @retval 0: Success.
int read_cafemol_dcd(std::ifstream&, std::vector<Conformation>&);

}
#endif

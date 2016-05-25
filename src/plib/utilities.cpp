/*!
************************************************************
@file utilities.cpp
@brief Defines lots of useful functions.

Definitions of useful functions.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-25 23:59
@copyright GNU Public License V3.0
************************************************************
*/

#include "utilities.hpp"

namespace pinang {

std::vector<std::string>& split_str(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split_str(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split_str(s, delim, elems);
  return elems;
}

}  // pinang

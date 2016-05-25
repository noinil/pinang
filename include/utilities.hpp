/*!
@file utilities.hpp
@brief Basic useful functions.

In this file a few useful handy functions are provided.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-26 00:01
@copyright GNU Public License V3.0
*/

#ifndef PINANG_UTILITIES_
#define PINANG_UTILITIES_

#include <string>
#include <sstream>
#include <vector>

namespace pinang {

//! @brief Split a long string into short strings (in a vector) by delim.
//! @param Original string.
//! @param Delim (separator).
//! @param Vector which is to store substrings.
//! @return Reference of vector which stores substrings.
std::vector<std::string>& split_str(const std::string&, char, std::vector<std::string>&);

//! @brief Split a long string into short strings (in a vector) by delim.
//! @param Original string.
//! @param Delim (separator).
//! @return Vector which stores substrings.
std::vector<std::string> split_str(const std::string&, char);

}

#endif

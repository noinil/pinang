/*!
************************************************************
@file selection.cpp
@brief Define functions of class Selection.

Definitions of member or friend functions of class Selection.

@author Cheng Tan (noinil@gmail.com)
@date 2016-05-25 22:58
@copyright GNU Public License V3.0
************************************************************
*/

#include <iostream>
#include <fstream>
#include "selection.hpp"
#include "utilities.hpp"

namespace pinang {

Selection::Selection()
{
  n_atom_ = 0;
  v_serial_.clear();
}

Selection::Selection(std::vector<int> v)
{
  v_serial_ = v;
  n_atom_ = v_serial_.size();
}

Selection::Selection(std::string inp_file_name, std::string keyword)
{
  v_serial_.clear();
  std::ifstream inp_file(inp_file_name.c_str());
  std::string inp_line;
  int check_flag = 1;
  std::string::size_type m;
  std::string::size_type n;
  std::string::size_type q;

  while (inp_file.good()) {
    std::getline(inp_file, inp_line);
    if (inp_file.fail())
      break;

    m = inp_line.find(keyword);
    if (m != 0) {
      continue;
    } else {
      check_flag = 0;
      n = inp_line.find(":");
      std::string select_str = inp_line.substr(n + 1);
      std::vector<std::string> sels = split_str(select_str, ',');
      int tmp_i = 0;
      int tmp_j = 0;
      for (const std::string& sel_word : sels){
        q = sel_word.find("to");
        if (q != std::string::npos){
          std::string sub1 = sel_word.substr(0, q);
          std::string sub2 = sel_word.substr(q + 2);
          tmp_i = std::stoi(sub1);
          tmp_j = std::stoi(sub2);
          for (int j = tmp_i; j <= tmp_j; j++) {
            v_serial_.push_back(j - 1);
          }
        } else {
          tmp_i = std::stoi(sel_word);
          v_serial_.push_back(tmp_i-1);
        }
      }
    }
  }
  inp_file.close();

  if (check_flag) {
    std::cout << "ERROR! Keyword " << keyword << " not found! --- in file selection.cpp. \n";
    exit(EXIT_SUCCESS);
  }

  n_atom_ = v_serial_.size();

}

void Selection::reset()
{
  n_atom_ = 0;
  v_serial_.clear();
}

int Selection::set_selection(std::vector<int> v)
{
  n_atom_ = v.size();
  v_serial_ = v;
  return 1;
}

int Selection::get_selection(int n)
{
  if (n >= n_atom_ || n < 0)
  {
    std::cout << " ~            PINANG :: selection.hpp       ~ " << "\n";
    std::cerr << " ERROR: Atom index out of range in Selection. " << "\n";
    exit(EXIT_SUCCESS);
  } else {
    return v_serial_[n];
  }
}


int Selection::set_selection(std::string inp_file_name, std::string keyword)
{
  v_serial_.clear();
  std::ifstream inp_file(inp_file_name.c_str());
  std::string inp_line;
  int check_flag = 1;
  std::string::size_type m;
  std::string::size_type n;
  std::string::size_type q;

  while (inp_file.good()) {
    std::getline(inp_file, inp_line);
    if (inp_file.fail())
      break;

    m = inp_line.find(keyword);
    if (m != 0) {
      continue;
    } else {
      check_flag = 0;
      n = inp_line.find(":");
      std::string select_str = inp_line.substr(n + 1);
      std::vector<std::string> sels = split_str(select_str, ',');
      int tmp_i = 0;
      int tmp_j = 0;
      for (const std::string& sel_word : sels){
        q = sel_word.find("to");
        if (q != std::string::npos){
          std::string sub1 = sel_word.substr(0, q);
          std::string sub2 = sel_word.substr(q + 2);
          tmp_i = std::stoi(sub1);
          tmp_j = std::stoi(sub2);
          for (int j = tmp_i; j <= tmp_j; j++) {
            v_serial_.push_back(j - 1);
          }
        } else {
          tmp_i = std::stoi(sel_word);
          v_serial_.push_back(tmp_i-1);
        }
      }
    }
  }
  inp_file.close();

  if (check_flag) {
    std::cout << "ERROR! Keyword " << keyword << " not found! --- in file selection.cpp. \n";
    exit(EXIT_SUCCESS);
  }

  n_atom_ = v_serial_.size();
  return 0;
}

}  // pinang

/*!
  @file selection.hpp
  @brief Definition of class Selection.

  In this file class Selection is defined.  Each selection contains a list of
  atom/residue indeces.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-05-25 22:37
  @copyright GNU Public License V3.0
*/

#ifndef PINANG_SELECTION_H
#define PINANG_SELECTION_H

#include <vector>
#include <string>

namespace pinang {

/*!
  @brief A set of index of selected atoms/residues.

  The class Selection consists of a set of indeces, which are of atoms/residues
  selected from a molecule.
*/
class Selection
{
 public:
  //! @brief Create an "empty" Selection object.
  //! @return A Selection object.
  Selection();
  //! @brief Create a Selection object based on a set of coordinates.
  //! @param A set of Vec3d coordinates.
  //! @return A Selection object.
  Selection(std::vector<int>);

  //! @brief Create a Selection object by reading KEYWORD from input file.
  //! @param An input file.
  //! @param KEYWORD.
  //! @return A Selection object.
  Selection(std::string, std::string);

  virtual ~Selection() {v_serial_.clear();}

  //! @brief Reset properties of Selection.
  void reset();

  //! @brief Get number of atoms/residues in Selection.
  //! @return Number of atoms/residues in Selection.
  int get_size() const { return n_atom_; }

  //! @brief Set Selection based on a set of coordinates.
  //! @param A set of serial numbers.
  //! @return Status of setting up the Selection.
  //! @retval 1: Failure.
  //! @retval 0: Success.
  int set_selection(std::vector<int>);

  //! @brief Set Selection by reading from input file.
  //! @param An input file.
  //! @param KEYWORD.
  //! @return Status of setting up the Selection.
  //! @retval 1: Failure.
  //! @retval 0: Success.
  int set_selection(std::string, std::string);

  //! @brief Get a index number of atom/residue from selection.
  //! @param Serial number of the selection.
  //! @return Index of atom/residue.
  int get_selection(int) const;

 protected:
  std::vector<int> v_serial_;       //!< A set of indecies (serial numbers).
  int n_atom_;                      //!< Number of coordinates in a selection.
};


}

#endif

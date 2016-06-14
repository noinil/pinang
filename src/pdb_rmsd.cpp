/*!
  @file pdb_rmsd.cpp
  @brief Calculate RMSD for structures from PDB.

  Read PDB, extract coordinates, calculate RMSD.

  @author Cheng Tan (noinil@gmail.com)
  @date 2016-06-14 15:45
  @copyright GNU Public License V3.0
*/

#include <fstream>
#include <unistd.h>
#include "conformation.hpp"
#include "selection.hpp"
#include "geometry.hpp"

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
  int opt;
  int in_flag = 0;
  std::string inpname = "selection.in";
  while ((opt = getopt(argc, argv, "i:h")) != -1) {
    switch (opt) {
      case 'i':
        inpname = optarg;
        in_flag = 1;
        break;
      case 'h':
        print_usage(argv[0]);
        break;
      default: /* '?' */
        print_usage(argv[0]);
    }
  }

  std::string crd_name1 = argv[argc - 2];
  std::string crd_name2 = argv[argc - 1];
  pinang::Conformation c1(crd_name1);
  pinang::Conformation c2(crd_name2);

  pinang::Selection ref1(inpname, "REFERENCE1");
  pinang::Selection ref2(inpname, "REFERENCE2");
  pinang::Selection cmp1(inpname, "CALCULATE1");
  pinang::Selection cmp2(inpname, "CALCULATE2");

  if (ref1.get_size() != ref2.get_size())
  {
    cout << " Error: Inconsistent Selection size. References groups. \n";
    print_usage(argv[0]);
  }
  if (cmp1.get_size() != cmp2.get_size())
  {
    cout << " Error: Inconsistent Selection size. RMSD Calculation groups. \n";
    print_usage(argv[0]);
  }

  pinang::Group ref_group1(c1, ref1);
  pinang::Group ref_group2(c2, ref2);
  pinang::Group cmp_group1(c1, cmp1);
  pinang::Group cmp_group2(c2, cmp2);

  pinang::Transform t;
  find_transform(ref_group1, ref_group2, t);
  pinang::Group trans_group1;
  trans_group1 = t.apply(cmp_group1);
  cout << get_rmsd(trans_group1, cmp_group2) << "\n";

  return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: " << s
            << "\n\t -i some.in struct1.crd struct2.crd\n\t [-h] \n";
  exit(EXIT_SUCCESS);
}

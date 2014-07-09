#include "read_dcd.h"

#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
    std::ifstream ifile(argv[1], std::ifstream::binary);
    std::vector<pinang::Conformation> conformations;

    std::ofstream ofile("test.coor");

    pinang::read_dcd(ifile, conformations);
    int tmp;
    int nframe = conformations.size();
    std::cout << nframe << std::endl;
    for (int i= 0; i < nframe; i++) {
        // std::cout << "Frame: " << i
        //           << " an: " << conformations[i].m_size()
        //           << std::endl;
        ofile << "Frame: " << i
                  << std::endl;
        for (int j = 0; j < conformations[i].m_size(); j++) {
            ofile << conformations[i].atom(j)
                  << std::endl;
        }
        // ofile << conformations[i].atom(1) << std::endl;
    }

    return 0;
}

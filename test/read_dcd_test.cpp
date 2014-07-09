#include "read_dcd.h"

#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
    std::ifstream ifile(argv[1]);
    std::vector<pinang::Vec3d> test;

    read_dcd(ifile, test);
    return 0;
}

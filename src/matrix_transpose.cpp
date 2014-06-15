#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

int main(int argc, char *argv[])
{
    std::vector < double > mat_line;
    std::vector < std::vector < double > > mat;

    std::cout << " ========================================================== "
              << std::endl;
    std::cout << " Usage: "
              << argv[0]
              << " -f some.dat [-o out.dat] [-h]" << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt, mod_index = 0;
    int mod_flag = 0;
    int in_flag = 0;
    std::string infilename = "some.dat";
    std::string outfilename = "out.dat";

    while ((opt = getopt(argc, argv, "o:f:h")) != -1) {
        switch (opt) {
        case 'o':
            outfilename = optarg;
            break;
        case 'f':
            infilename = optarg;
            in_flag = 1;
            break;
        case 'h':
            std::cout << " Help: Usage: "
                      << argv[0]
                      << " -f some.dat [-o out.dat] [-h]" << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dat [-o out.dat] [-h]" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " Usage: "
                  << argv[0]
                  << " -f some.dat [-o out.dat] [-h]" << std::endl;
        exit(EXIT_SUCCESS);
    }
    std::ifstream ifile(infilename.c_str());

    std::ofstream ofile(outfilename.c_str());

    double t;
    std::string file_line;
    std::istringstream tmp_sstr;

    if (!ifile.is_open())
    {
        std::cerr << " ERROR: Cannot read file: " << infilename << std::endl;
        exit(EXIT_SUCCESS);
    }

    while (ifile.good()) {
        std::getline(ifile, file_line);
        if (ifile.fail())
            break;
        tmp_sstr.str ( file_line);
        while (tmp_sstr.good()) {
            tmp_sstr >> t;
            if (tmp_sstr.fail())
                break;
            mat_line.push_back(t);
        }
        mat.push_back(mat_line);
        // TEST: test mat reading and vector push_back;
        // std::cout << mat.size() << std::endl;
        mat_line.clear();
        tmp_sstr.clear();
    }

    int i,j;
    for (i = 0; i < mat[0].size(); i++) {
        for (j = 0; j < mat.size(); j++) {
            std::cout << std::setw(10)
                      << std::setiosflags(std::ios_base::fixed)
                      << std::setprecision(2)
                      << mat[j][i];
            ofile << std::setw(10)
                  << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(2)
                  << mat[j][i];
        }
        std::cout << std::endl;
        ofile << std::endl;
    }

    ifile.close();
    ofile.close();

    return 0;
}

#include "read_dcd.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <boost/algorithm/string.hpp>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~           PINANG DCD distances calculation             ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt;
    int inp_flag = 0;
    int dcd_flag = 0;
    int top_flag = 0;

    std::string dcd_name = "some.dcd";
    std::string top_name = "some.top";
    std::string inp_name = "some.inp";
    std::string dis_name = "some.dis";

    while ((opt = getopt(argc, argv, "f:s:i:o:h")) != -1) {
        switch (opt) {
        case 'f':
            dcd_name = optarg;
            dcd_flag = 1;
            break;
        case 's':
            top_name = optarg;
            top_flag = 1;
            break;
        case 'i':
            inp_name = optarg;
            inp_flag = 1;
            break;
        case 'o':
            dis_name = optarg;
            break;
        case 'h':
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.dis] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.dis] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (dcd_flag == 0)
    {
        std::cout << " ERROR: Please provide the DCD file. " << std::endl;
        exit(EXIT_SUCCESS);
    }
    if (top_flag == 0)
    {
        std::cout << " ERROR: Please provide the top file. " << std::endl;
        exit(EXIT_SUCCESS);
    }
    if (inp_flag == 0)
    {
        std::cout << " ERROR: Please provide the input file. " << std::endl;
        exit(EXIT_SUCCESS);
    }
    // -------------------------------------------------------------------------

    std::ifstream dcd_file(dcd_name.c_str(), std::ifstream::binary);
    std::ifstream inp_file(inp_name.c_str());
    std::ifstream top_file(top_name.c_str());
    std::ofstream dis_file(dis_name.c_str());

    int N_atom = 0;
    std::string str0 = "particles";
    std::string inp_line;
    while (top_file.good()) {
        std::getline(top_file, inp_line);

        if (top_file.fail())
        {
            break;
        }

        std::size_t found = inp_line.find(str0);
        if (found!=std::string::npos){
            std::cout << " - TOP file: total number "
                      << inp_line << std::endl;
            std::string stmp;
            std::istringstream tmp_sstr;
            tmp_sstr.str ( inp_line );

            tmp_sstr >> stmp  >> stmp  >> stmp
                     >> N_atom;
            break;
        }
    }
    if (N_atom == 0)
    {
        std::cout << " ERROR: No particles found in top file. " << std::endl;
        exit(EXIT_SUCCESS);
    }


    // ==================== input particle index ====================
    std::istringstream tmp_sstr;
    int tmp_a, tmp_b;
    int flg_grp_1 = 0;
    int flg_grp_2 = 0;
    std::string tmp_str;

    while (inp_file.good()) {
        std::getline(inp_file, inp_line);

        if (inp_file.fail())
        {
            break;
        }

        tmp_str = inp_line.substr(0,7);

        if (tmp_str == "GROUP1:")
        {
            flg_grp_1 = 1;
            inp_line.erase(0,7);
            std::vector<std::string> strs;
            boost::split(strs, inp_line, boost::is_any_of(","));
            for (int i = 0; i < strs.size(); i++) {
                std::cout << strs[i] << std::endl;
            }
        }
    }
    if (flg_grp_1 == 0)
    {
        std::cout << " ERROR: GROUP 1 not found!" << std::endl;
    }
    inp_file.close();
    exit(EXIT_SUCCESS);

    // -------------------------------------------------------------------------
    // ---------- Reading DCD ----------
    std::vector<pinang::Conformation> conformations;

    pinang::read_dcd(dcd_file, conformations);
    int nframe = conformations.size();

    if (nframe == 0)
    {
        std::cout << " ERROR: Empty DCD file!  Please check! " << std::endl;
        return 1;
    }

    if (N_atom != conformations[0].m_size())
    {
        std::cout << " ERROR: Particle number don't match in top and dcd! "
                  << " Please check! " << std::endl;
        return 1;
    }

    /*                  _         _
    //  _ __ ___   __ _(_)_ __   | | ___   ___  _ __
    // | '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \
    // | | | | | | (_| | | | | | | | (_) | (_) | |_) |
    // |_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/
    //                                         |_|
    */

    dcd_file.close();
    top_file.close();
    dis_file.close();

    return 0;
}

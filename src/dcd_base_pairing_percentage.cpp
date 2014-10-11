#include "read_dcd.h"
#include "topology.h"

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
    std::cout << " ~        PINANG DCD DNA base pairing calculation         ~ "
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
    std::string bpn_name = "some.bpn";

    while ((opt = getopt(argc, argv, "f:s:i:o:c:h")) != -1) {
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
            bpn_name = optarg;
            break;
        case 'c':
            pinang::g_cutoff = atof(optarg);
            break;
        case 'h':
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.bpn] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.bpn] [-h]"
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
    std::ofstream bpn_file(bpn_name.c_str());

    pinang::Topology top(top_name);
    if (top.m_size() == 0)
    {
        std::cout << " ERROR: No particles found in top file. " << std::endl;
        exit(EXIT_SUCCESS);
    }

    // ==================== input particle index ====================
    std::istringstream tmp_sstr;
    int flg_grp_1 = 0;
    int flg_grp_2 = 0;
    std::string tmp_str;
    std::string inp_line;

    std::vector<int> atom_group1_idx;
    std::vector<int> atom_group2_idx;

    while (inp_file.good()) {
        std::getline(inp_file, inp_line);

        if (inp_file.fail())
        {
            break;
        }

        tmp_str = inp_line.substr(0,7);

        if (tmp_str == "GROUP1:")
        {
            std::string tmp_s;
            std::istringstream tmp_sstr;
            flg_grp_1 = 1;
            inp_line.erase(0,7);
            std::vector<std::string> strs;
            boost::split(strs, inp_line, boost::is_any_of(","));
            for (int i = 0; i < strs.size(); i++) {
                tmp_s = strs[i];
                std::size_t found = tmp_s.find("to");
                std::size_t found1 = tmp_s.find("every");
                if (found!=std::string::npos && found1 == std::string::npos){
                    int tmp_i = 0;
                    int tmp_j = 0;
                    tmp_s.erase(found, 2);
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    tmp_sstr >> tmp_j;
                    for (int j = tmp_i; j <= tmp_j; j++) {
                        atom_group1_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else if (found!=std::string::npos && found1!=std::string::npos){
                    int tmp_i = 0;
                    int tmp_j = 0;
                    int tmp_k = 0;
                    tmp_s.erase(found, 2);
                    tmp_s.erase(found1-2, 5);
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    tmp_sstr >> tmp_j;
                    tmp_sstr >> tmp_k;
                    if (tmp_k > 0) {
                        for (int j = tmp_i; j <= tmp_j; j+=tmp_k)
                            atom_group1_idx.push_back(j-1);
                    }
                    if (tmp_k < 0) {
                        for (int j = tmp_i; j >= tmp_j; j+=tmp_k)
                            atom_group1_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else {
                    int tmp_i = 0;
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    atom_group1_idx.push_back(tmp_i-1);
                    tmp_sstr.clear();
                }
            }
        }
        if (tmp_str == "GROUP2:")
        {
            std::string tmp_s;
            std::istringstream tmp_sstr;
            flg_grp_2 = 1;
            inp_line.erase(0,7);
            std::vector<std::string> strs;
            boost::split(strs, inp_line, boost::is_any_of(","));
            for (int i = 0; i < strs.size(); i++) {
                tmp_s = strs[i];
                std::size_t found = tmp_s.find("to");
                std::size_t found1 = tmp_s.find("every");
                if (found!=std::string::npos && found1 == std::string::npos){
                    int tmp_i = 0;
                    int tmp_j = 0;
                    tmp_s.erase(found, 2);
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    tmp_sstr >> tmp_j;
                    for (int j = tmp_i; j <= tmp_j; j++) {
                        atom_group2_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else if (found!=std::string::npos && found1 != std::string::npos){
                    int tmp_i = 0;
                    int tmp_j = 0;
                    int tmp_k = 0;
                    tmp_s.erase(found, 2);
                    tmp_s.erase(found1-2, 5);
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    tmp_sstr >> tmp_j;
                    tmp_sstr >> tmp_k;
                    if (tmp_k > 0)
                    {
                        for (int j = tmp_i; j <= tmp_j; j+=tmp_k)
                            atom_group2_idx.push_back(j-1);
                    }
                    if (tmp_k < 0)
                    {
                        for (int j = tmp_i; j >= tmp_j; j+=tmp_k)
                            atom_group2_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else {
                    int tmp_i = 0;
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    atom_group2_idx.push_back(tmp_i-1);
                    tmp_sstr.clear();
                }
            }
        }
    }
    // for (int i = 0; i < atom_group1_idx.size(); i++) {
    //     std::cout << atom_group1_idx[i] << "   ";
    //     std::cout << atom_group2_idx[i] << std::endl;
    // }

    if (flg_grp_1 == 0)
    {
        std::cout << " ERROR: GROUP 1 not found!" << std::endl;
        return 1;
    }
    if (flg_grp_2 == 0)
    {
        std::cout << " ERROR: GROUP 2 not found!" << std::endl;
        return 1;
    }
    if (atom_group1_idx.size() != atom_group2_idx.size())
    {
        std::cout << " ERROR: GROUP 1 and 2 not corresponding!" << std::endl;
        return 1;
    }
    inp_file.close();

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

    if (top.m_size() != conformations[0].m_size())
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
    pinang::Vec3d com1, com2;
    double dist = 0;
    for (int i= 0; i < nframe; i++) {
        pinang::Vec3d v(0,0,0);
        int bp_number = 0;
        int k = 0, l = 0;

        for (int j = 0; j < atom_group1_idx.size(); j++) {
            k = atom_group1_idx[j];
            l = atom_group2_idx[j];
            v = conformations[i].atom(k) - conformations[i].atom(l);
            dist = v.norm();
            if (dist <= pinang::g_cutoff)
                bp_number++;
        }

        bpn_file << std::setw(6) << i
                 << "   " << std::setw(8) << bp_number
                 << std::endl; // Output the distance!
    }

    dcd_file.close();
    bpn_file.close();

    return 0;
}
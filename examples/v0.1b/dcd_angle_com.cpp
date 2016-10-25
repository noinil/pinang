#include "read_cafemol_dcd.hpp"
#include "topology.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <boost/algorithm/string.hpp>

using namespace std;

void print_usage(char* s);

int main(int argc, char *argv[])
{
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
          print_usage(argv[0]);
           break;
        default: /* '?' */
          print_usage(argv[0]);
       }
    }

    if (dcd_flag == 0)
    {
        std::cout << " ERROR: Please provide the DCD file. " << "\n";
        exit(EXIT_SUCCESS);
    }
    if (top_flag == 0)
    {
        std::cout << " ERROR: Please provide the top file. " << "\n";
        exit(EXIT_SUCCESS);
    }
    if (inp_flag == 0)
    {
        std::cout << " ERROR: Please provide the input file. " << "\n";
        exit(EXIT_SUCCESS);
    }
    // -------------------------------------------------------------------------

    std::ifstream dcd_file(dcd_name.c_str(), std::ifstream::binary);
    std::ifstream inp_file(inp_name.c_str());
    std::ofstream dis_file(dis_name.c_str());

    pinang::Topology top(top_name);
    if (top.get_size() == 0)
    {
        std::cout << " ERROR: No particles found in top file. " << "\n";
        exit(EXIT_SUCCESS);
    }

    // ==================== input particle index ====================
    std::istringstream tmp_sstr;
    int flg_grp_1 = 0;
    int flg_grp_2 = 0;
    int flg_grp_3 = 0;
    std::string tmp_str;
    std::string inp_line;

    std::vector<int> atom_group1_idx;
    std::vector<int> atom_group2_idx;
    std::vector<int> atom_group3_idx;

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
            for (int i = 0; i < int(strs.size()); i++) {
                tmp_s = strs[i];
                std::size_t found = tmp_s.find("to");
                if (found!=std::string::npos){
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
            for (int i = 0; i < int(strs.size()); i++) {
                tmp_s = strs[i];
                std::size_t found = tmp_s.find("to");
                if (found!=std::string::npos){
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
                } else {
                    int tmp_i = 0;
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    atom_group2_idx.push_back(tmp_i-1);
                    tmp_sstr.clear();
                }
            }
        }
        if (tmp_str == "GROUP3:")
        {
            std::string tmp_s;
            std::istringstream tmp_sstr;
            flg_grp_3 = 1;
            inp_line.erase(0,7);
            std::vector<std::string> strs;
            boost::split(strs, inp_line, boost::is_any_of(","));
            for (int i = 0; i < int(strs.size()); i++) {
                tmp_s = strs[i];
                std::size_t found = tmp_s.find("to");
                if (found!=std::string::npos){
                    int tmp_i = 0;
                    int tmp_j = 0;
                    tmp_s.erase(found, 2);
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    tmp_sstr >> tmp_j;
                    for (int j = tmp_i; j <= tmp_j; j++) {
                        atom_group3_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else {
                    int tmp_i = 0;
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    atom_group3_idx.push_back(tmp_i-1);
                    tmp_sstr.clear();
                }
            }
        }

    }
    if (flg_grp_1 == 0)
    {
        std::cout << " ERROR: GROUP 1 not found!" << "\n";
    }
    if (flg_grp_2 == 0)
    {
        std::cout << " ERROR: GROUP 2 not found!" << "\n";
    }
    if (flg_grp_3 == 0)
    {
        std::cout << " ERROR: GROUP 2 not found!" << "\n";
    }
    inp_file.close();

    // -------------------------------------------------------------------------
    // ---------- Reading DCD ----------
    std::vector<pinang::Conformation> conformations;

    pinang::read_cafemol_dcd(dcd_file, conformations);
    int nframe = conformations.size();

    if (nframe == 0)
    {
        std::cout << " ERROR: Empty DCD file!  Please check! " << "\n";
        return 1;
    }

    if (top.get_size() != conformations[0].get_size())
    {
        std::cout << " ERROR: Particle number don't match in top and dcd! "
                  << " Please check! " << "\n";
        return 1;
    }

    /*                  _         _
    //  _ __ ___   __ _(_)_ __   | | ___   ___  _ __
    // | '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \
    // | | | | | | (_| | | | | | | | (_) | (_) | |_) |
    // |_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/
    //                                         |_|
    */
    pinang::Vec3d com1, com2, com3;
    pinang::Vec3d dnav1, dnav2;
    double angle = 0;
    for (int i= 0; i < nframe; i++) {
        pinang::Vec3d v(0,0,0);
        int k = 0;
        double total_mass = 0;

        for (int j = 0; j < int(atom_group1_idx.size()); j++) {
            k = atom_group1_idx[j];
            v = v + conformations[i].get_coordinate(k) * top.get_particle(k).get_mass();
            total_mass += top.get_particle(k).get_mass();
        }
        com1 = v * (1/total_mass);

        v = v * 0;
        total_mass = 0;
        for (int j = 0; j < int(atom_group2_idx.size()); j++) {
            k = atom_group2_idx[j];
            v = v + conformations[i].get_coordinate(k) * top.get_particle(k).get_mass();
            total_mass += top.get_particle(k).get_mass();
        }
        com2 = v * (1/total_mass);

        v = v * 0;
        total_mass = 0;
        for (int j = 0; j < int(atom_group3_idx.size()); j++) {
            k = atom_group3_idx[j];
            v = v + conformations[i].get_coordinate(k) * top.get_particle(k).get_mass();
            total_mass += top.get_particle(k).get_mass();
        }
        com3 = v * (1/total_mass);

        dnav1 = com2 - com1;
        dnav2 = com3 - com2;
        angle = vec_angle_deg(dnav1, dnav2);

        dis_file << std::setw(6) << i
                 << "   " << std::setw(8) << angle
                 << "\n"; // Output the distance!
    }

    dcd_file.close();
    dis_file.close();

    return 0;
}

void print_usage(char* s)
{
  std::cout << " Usage: "
            << s
            << " -f some.dcd -s some.top -i some.inp [-o some.dis] [-h]"
            << "\n";
  exit(EXIT_SUCCESS);
}

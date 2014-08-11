#include "read_dcd.h"
#include "topology.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <numeric>
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~      PINANG DCD persistence length calculation         ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt;
    int inp_flag = 0;
    int dcd_flag = 0;
    int top_flag = 0;
    int beg_flag = 0;
    int end_flag = 0;
    int begin_frame = 0;
    int end_frame = 0;

    std::string dcd_name = "some.dcd";
    std::string top_name = "some.top";
    std::string inp_name = "some.inp";
    std::string lp_name = "some.lp";

    while ((opt = getopt(argc, argv, "f:s:i:o:b:e:h")) != -1) {
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
            lp_name = optarg;
            break;
        case 'b':
            begin_frame = atoi(optarg);
            beg_flag = 1;
            break;
        case 'e':
            end_frame = atoi(optarg);
            end_flag = 1;
            break;
        case 'h':
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.lp]\n"
                      << " [-b] beginning_frame [-e] ending_frame [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.lp]\n"
                      << " [-b] beginning_frame [-e] ending_frame [-h]"
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

    std::ifstream dcd_file(dcd_name.c_str(), std::ifstream::binary);
    std::ifstream inp_file(inp_name.c_str());
    std::ofstream lp_file(lp_name.c_str());

    // ==================== topology file read in ====================
    pinang::Topology top(top_name);
    std::cout << " (i.e. " << (top.m_size()+2)/6 << " bp.)" << std::endl;
    if (top.m_size() == 0)
    {
        std::cout << " ERROR: No particles found in top file. " << std::endl;
        exit(EXIT_SUCCESS);
    }

    // -------------------------------------------------------------------------

    // ==================== input particle index ====================
    std::istringstream tmp_sstr;
    int flg_grp_1 = 0;
    int flg_grp_2 = 0;
    std::string tmp_str;
    std::string inp_line;

    std::vector<int> phosphate_strand1_idx;
    std::vector<int> phosphate_strand2_idx;
    std::vector<int> base_strand1_idx;
    std::vector<int> base_strand2_idx;

    while (inp_file.good()) {
        std::getline(inp_file, inp_line);

        if (inp_file.fail())
        {
            break;
        }

        tmp_str = inp_line.substr(0,8);

        if (tmp_str == "STRAND1:")
        {
            std::string tmp_s;
            std::istringstream tmp_sstr;
            flg_grp_1 = 1;
            inp_line.erase(0,8);
            std::vector<std::string> strs;
            boost::split(strs, inp_line, boost::is_any_of(","));
            for (int i = 0; i < strs.size(); i++) {
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
                        if (top.particle(j-1).atom_name()=="P  ")
                            phosphate_strand1_idx.push_back(j-1);
                        if (top.particle(j-1).atom_name()=="B  ")
                            base_strand1_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else {
                    int tmp_i = 0;
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    if (top.particle(tmp_i-1).atom_name()=="P  ")
                        phosphate_strand1_idx.push_back(tmp_i-1);
                    if (top.particle(tmp_i-1).atom_name()=="B  ")
                        base_strand1_idx.push_back(tmp_i-1);
                    tmp_sstr.clear();
                }
            }
        }

        if (tmp_str == "STRAND2:")
        {
            std::string tmp_s;
            std::istringstream tmp_sstr;
            flg_grp_2 = 1;
            inp_line.erase(0,8);
            std::vector<std::string> strs;
            boost::split(strs, inp_line, boost::is_any_of(","));
            for (int i = 0; i < strs.size(); i++) {
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
                        if (top.particle(j-1).atom_name()=="P  ")
                            phosphate_strand2_idx.push_back(j-1);
                        if (top.particle(j-1).atom_name()=="B  ")
                            base_strand2_idx.push_back(j-1);
                    }
                    tmp_sstr.clear();
                } else {
                    int tmp_i = 0;
                    tmp_sstr.str(tmp_s);
                    tmp_sstr >> tmp_i;
                    if (top.particle(tmp_i-1).atom_name()=="P  ")
                        phosphate_strand2_idx.push_back(tmp_i-1);
                    if (top.particle(tmp_i-1).atom_name()=="B  ")
                        base_strand2_idx.push_back(tmp_i-1);
                    tmp_sstr.clear();
                }
            }
        }
    }
    if (flg_grp_1 == 0)
    {
        std::cout << " ERROR: STRAND 1 not found!" << std::endl;
    }
    if (flg_grp_2 == 0)
    {
        std::cout << " ERROR: STRAND 2 not found!" << std::endl;
    }
    inp_file.close();
    if (phosphate_strand1_idx.size() != phosphate_strand2_idx.size())
    {
        std::cout << " ERROR: STRAND 1 Phosphate != STRAND 2 Phosphate!"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    if (base_strand1_idx.size() != base_strand2_idx.size())
    {
        std::cout << " ERROR: STRAND 1 and STRAND 2 bases do not match!"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    if (base_strand1_idx.size() != phosphate_strand1_idx.size() + 1)
    {
        std::cout << " ERROR: Number of Bases should be \n"
                  << " one more than number of phosphates!"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }

    // -------------------------------------------------------------------------
    // ---------- Reading DCD ----------
    std::vector<pinang::Conformation> conformations;

    pinang::read_dcd(dcd_file, conformations);
    int nframe = conformations.size();
    std::cout << " Total " << nframe << " frames in dcd file." << std::endl;

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
    if (beg_flag == 0)
    {
        begin_frame = 0;
    }
    if (end_flag == 0)
    {
        end_frame = nframe;
    } else if (end_frame > nframe)
    {
        end_frame = nframe;
    }
    std::cout << " Calculating from " << begin_frame << "th frame to "
              << end_frame << "th frame..." << std::endl;
    /*                  _         _
    //  _ __ ___   __ _(_)_ __   | | ___   ___  _ __
    // | '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \
    // | | | | | | (_| | | | | | | | (_) | (_) | |_) |
    // |_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/
    //                                         |_|
    */
    const double p_D_phosphate = 17.45; // diameter of helix cross section;

    // Persistence_length Definition 1: lp = <Re^2>/2L;
    double persistence_length_1 = 0;
    // Persistence_length Definition 2: <u(0) * u(s)> = e^(-s/lp);
    double persistence_length_2 = 0;
    // Contour length: \sum ( base_rise );
    double contour_length = 0;  // This is an average of lc<frames>;

    std::vector<double> R_e_sqr; // Re^2 for each frame;
    std::vector<double> lc;      // contour length for each frame;
    std::vector<double> u_u_50;  // u(0) * u(50%);
    std::vector<double> u_u_60;  // u(0) * u(60%);
    std::vector<double> u_u_70;  // u(0) * u(70%);
    std::vector<double> u_u_80;  // u(0) * u(80%);
    std::vector<double> u_u_90;  // u(0) * u(90%);
    std::vector<double> u_u_lc;  // u(0) * u(100%);

    for (int h = begin_frame; h < end_frame; h++) {
        std::vector<pinang::Vec3d> curve1_nodes; // "P" atoms in the 1st backbone;
        std::vector<pinang::Vec3d> curve2_nodes; // "P " atoms in the 2nd backbone;
        std::vector<pinang::Vec3d> cross_vecs;
        std::vector<pinang::Vec3d> base_positions;
        std::vector<pinang::Vec3d> axis_nodes;
        std::vector<pinang::Vec3d> axis_directions; // Axis direction vectors;
        std::vector<double> base_rise;
        pinang::Conformation confm = conformations[h];

        double d_ee = 0;        // end-end distance ^ 2;
        double d_lc = 0;        // sum of base rise;

        int len = phosphate_strand1_idx.size();
        for (int i = 0; i < len; i++) {
            int j = phosphate_strand1_idx[i];
            curve1_nodes.push_back(confm.atom(j));

            int k = phosphate_strand2_idx[i];
            curve2_nodes.push_back(confm.atom(k));
        }
        for (int i = 0; i < len + 1; i++) {
            int j = base_strand1_idx[i];
            base_positions.push_back(confm.atom(j));
        }

        pinang::Vec3d t_tmp(0,0,0);
        for (int i = 0; i < len; i++) {
            t_tmp = curve2_nodes[len-1-i] - curve1_nodes[i]; // for cross_vectors;
            cross_vecs.push_back(t_tmp);
        }

        pinang::Vec3d H(0,0,0);     // for use in the algorithm of axis calculation;
        pinang::Vec3d P1(0,0,0);    // for use in the algorithm of axis calculation;
        pinang::Vec3d P2(0,0,0);    // for use in the algorithm of axis calculation;
        pinang::Vec3d v1(0,0,0);    // for use in the algorithm of axis calculation;
        pinang::Vec3d v2(0,0,0);    // for use in the algorithm of axis calculation;
        pinang::Vec3d n1(0,0,0);    // for use in the algorithm of axis calculation;

        double d = 0;               // base rise tmp value;
        for (int i = 0; i < len; i++) {
            if (i == len-1) {
                v1 = cross_vecs[i-1];
                v2 = cross_vecs[i];
                t_tmp = v1^v2;
                H = t_tmp * (1/t_tmp.norm()); // H is important, the local direction;
                P1 = curve1_nodes[i] + (v2 * 0.5); // midpoint of cross_vecs;
                t_tmp = H ^ v2;
                n1 = t_tmp * (1/t_tmp.norm());
            } else {
                v1 = cross_vecs[i];
                v2 = cross_vecs[i+1];
                t_tmp = v1^v2;
                H = t_tmp * (1/t_tmp.norm()); // H is important, the local direction;
                P1 = curve1_nodes[i] + (v1 * 0.5); // midpoint of cross_vecs;
                t_tmp = H ^ v1;
                n1 = t_tmp * (1/t_tmp.norm());
            }
            axis_directions.push_back(H); // axis direction vectors!

            // r = axis center deviation!  This is just a empirical estimation!
            double r = 2.12 * v1.norm() / p_D_phosphate;
            t_tmp = P1 + (n1 * r);
            axis_nodes.push_back(t_tmp);

            // Calculation of base rise!
            P1 = base_positions[i];
            P2 = base_positions[i+1];
            d = (P2 - P1) * H;  // base rise
            base_rise.push_back(d);
        }

        int last = axis_nodes.size()-1;
        n1 = axis_nodes[0] - axis_nodes[last];
        d_ee = n1.sqr_norm();
        R_e_sqr.push_back(d_ee);

        for (int j = 0; j < base_rise.size(); j++) {
            d_lc += base_rise[j];
        }
        lc.push_back(d_lc);

        double u_u_tmp = 0;
        int s = int(len * 0.5);
        u_u_tmp = axis_directions[0] * axis_directions[s-1];
        u_u_50.push_back(u_u_tmp);
        s = int(len * 0.6);
        u_u_tmp = axis_directions[0] * axis_directions[s-1];
        u_u_60.push_back(u_u_tmp);
        s = int(len * 0.7);
        u_u_tmp = axis_directions[0] * axis_directions[s-1];
        u_u_70.push_back(u_u_tmp);
        s = int(len * 0.8);
        u_u_tmp = axis_directions[0] * axis_directions[s-1];
        u_u_80.push_back(u_u_tmp);
        s = int(len * 0.9);
        u_u_tmp = axis_directions[0] * axis_directions[s-1];
        u_u_90.push_back(u_u_tmp);

        u_u_tmp = axis_directions[0] * axis_directions[len-1];
        u_u_lc.push_back(u_u_tmp);

        // lp_file << std::setw(6) << i
        //         << "   " << std::setw(8) << persistence_length_2
        //         << std::endl;
    }

    // ==================== Calculating Def1 of Lp ====================
    // ---------- calculating contour length ----------
    double sum = std::accumulate(lc.begin(), lc.end(), 0.0);
    contour_length = sum / lc.size() + 3.43;
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " Contour length = "
              << std::setw(8) << contour_length
              << std::endl;

    // ---------- calculating <Re^2> ----------
    sum = std::accumulate(R_e_sqr.begin(), R_e_sqr.end(), 0.0);
    double mean = sum / R_e_sqr.size();
    persistence_length_1 = mean / (2 * contour_length);

    // ~~~~~~~~~~ output ~~~~~~~~~~
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " Def 1 of persistence length: "
              << " lp = <Re^2>/2L" << std::endl
              << " Result: " << std::setw(8) << persistence_length_1
              << std::endl;
    lp_file << " Def 1 of persistence length: "
            << " lp = <Re^2>/2L" << std::endl
            << " Result: " << std::setw(8) << persistence_length_1
            << std::endl;
    lp_file << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << std::endl;

    // ==================== Calculating Def2 of Lp ====================
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " Def 2 of persistence length: "
              << " <u(0) * u(s)> = e^(-s/lp)" << std::endl;

    lp_file <<  " Def 2 of persistence length: "
            << " <u(0) * u(s)> = e^(-s/lp)" << std::endl;
    lp_file << std::setw(8) << "s" << "  "
            << std::setw(8) << "-Ln<u*u>"
            << std::setw(8) << "l_p" << std::endl;

    double p2_5, p2_6, p2_7, p2_8, p2_9, p2_10; // all calculated persistence lengthes;
    sum = std::accumulate(u_u_50.begin(), u_u_50.end(), 0.0);
    mean = sum / u_u_50.size();
    p2_5 = - 0.5 * contour_length / log(mean);
    lp_file << std::setw(8) << contour_length * 0.5 << "  "
            << std::setw(8) << -log(mean)
            << std::setw(8) << p2_5
            << std::endl;
    std::cout << " Result (s=50% contour length): " << std::setw(8)
              << p2_5
              << std::endl;

    sum = std::accumulate(u_u_60.begin(), u_u_60.end(), 0.0);
    mean = sum / u_u_60.size();
    p2_6 = - 0.6 * contour_length / log(mean);
    lp_file << std::setw(8) << contour_length * 0.6 << "  "
            << std::setw(8) << -log(mean)
            << std::setw(8) << p2_6
            << std::endl;
    std::cout << " Result (s=60% contour length): " << std::setw(8)
              << p2_6
              << std::endl;

    sum = std::accumulate(u_u_70.begin(), u_u_70.end(), 0.0);
    mean = sum / u_u_70.size();
    p2_7 = - 0.7 * contour_length / log(mean);
    lp_file << std::setw(8) << contour_length * 0.7 << "  "
            << std::setw(8) << -log(mean)
            << std::setw(8) << p2_7
            << std::endl;
    std::cout << " Result (s=70% contour length): " << std::setw(8)
              << p2_7
              << std::endl;

    sum = std::accumulate(u_u_80.begin(), u_u_80.end(), 0.0);
    mean = sum / u_u_80.size();
    p2_8 = - 0.8*contour_length / log(mean);
    lp_file << std::setw(8) << contour_length * 0.8 << "  "
            << std::setw(8) << -log(mean)
            << std::setw(8) << p2_8
            << std::endl;
    std::cout << " Result (s=80% contour length): " << std::setw(8)
              << p2_8
              << std::endl;

    sum = std::accumulate(u_u_90.begin(), u_u_90.end(), 0.0);
    mean = sum / u_u_90.size();
    p2_9 = - 0.9*contour_length / log(mean);
    lp_file << std::setw(8) << contour_length * 0.9 << "  "
            << std::setw(8) << -log(mean)
            << std::setw(8) << p2_9
            << std::endl;
    std::cout << " Result (s=90% contour length): " << std::setw(8)
              << p2_9
              << std::endl;

    sum = std::accumulate(u_u_lc.begin(), u_u_lc.end(), 0.0);
    mean = sum / u_u_lc.size();
    p2_10 = - contour_length / log(mean);
    lp_file << std::setw(8) << contour_length << "  "
            << std::setw(8) << -log(mean)
            << std::setw(8) << p2_10
            << std::endl;
    std::cout << " Result (s=100% contour length): " << std::setw(8)
              << p2_10
              << std::endl;

    persistence_length_2 = (p2_5 + p2_6 + p2_7 + p2_8 + p2_9 + p2_10)/6;
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    lp_file << " Averaged from above:"
            << std::setw(8) << persistence_length_2
            << std::endl;
    std::cout << " Averaged from above:"
              << std::setw(8) << persistence_length_2
              << std::endl;

    // ending ------------------------------
    dcd_file.close();
    lp_file.close();

    return 0;
}

#include "read_cafemol_dcd.h"
#include "topology.h"
#include "constants.h"

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

void gen_spline_fit(const std::vector<pinang::Vec3d>&, int N,
                    std::vector<pinang::Vec3d>&, std::vector<pinang::Vec3d>&);

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~           PINANG DCD DNA bending curvature             ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt;
    int dt = 100;
    int inp_flag = 0;
    int dcd_flag = 0;
    int top_flag = 0;
    int beg_flag = 0;
    int end_flag = 0;
    int log_flag = 0;
    int begin_frame = 0;
    int end_frame = 0;

    std::string dcd_name = "some.dcd";
    std::string top_name = "some.top";
    std::string inp_name = "some.inp";
    std::string crv_name = "some.curve";
    std::string log_name = "curve.log";

    while ((opt = getopt(argc, argv, "f:s:i:o:l:b:d:e:h")) != -1) {
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
            crv_name = optarg;
            break;
        case 'd':
            dt = std::max(atoi(optarg), 1);
            break;
        case 'l':
            log_name = optarg;
            log_flag = 1;
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
                      << " -f some.dcd -s some.top -i some.inp [-o some.curve]\n"
                      << " [-b] beginning_frame [-e] ending_frame [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.dcd -s some.top -i some.inp [-o some.curve]\n"
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
    std::ofstream crv_file(crv_name.c_str());
    std::ofstream log_file(log_name.c_str());

    // ==================== topology file read in ====================
    pinang::Topology top(top_name);
    if (top.get_size() == 0)
    {
        std::cout << " ERROR: No particles found in top file. " << std::endl;
        exit(EXIT_SUCCESS);
    }

    // -------------------------------------------------------------------------
    // ---------- Reading DCD ----------
    std::vector<pinang::Conformation> conformations;

    pinang::read_cafemol_dcd(dcd_file, conformations);
    int nframe = conformations.size();
    std::cout << " Total " << nframe << " frames in dcd file." << std::endl;

    if (nframe == 0)
    {
        std::cout << " ERROR: Empty DCD file!  Please check! " << std::endl;
        return 1;
    }

    if (top.get_size() != conformations[0].get_size())
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



    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // main
    //  _ _                 _      _ _
    // ( | )_ __ ___   __ _(_)_ __( | )
    //  V V| '_ ` _ \ / _` | | '_ \V V
    //     | | | | | | (_| | | | | |
    //     |_| |_| |_|\__,_|_|_| |_|
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    std::vector < std::vector <double> > curvature_all;

    /*      _                   _
    //     (_)_ __  _ __  _   _| |_
    //     | | '_ \| '_ \| | | | __|
    //     | | | | | |_) | |_| | |_
    //     |_|_| |_| .__/ \__,_|\__|
    //              |_|
    */
    int i = 0;
    int j = 0;
    double angle_lim = pinang::k_pi / 3;   // 60 degree;
    double angle_lim2 = pinang::k_pi / 15;   // 12 degree;
    double pi_over_36 = pinang::k_pi / 36; // 5 degree;
    double pi_over_60 = pinang::k_pi / 60; // 3 degree;

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

    /*                  _         _
    //  _ __ ___   __ _(_)_ __   | | ___   ___  _ __
    // | '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \
    // | | | | | | (_| | | | | | | | (_) | (_) | |_) |
    // |_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/
    //                                         |_|
    */
    int step = (end_frame - begin_frame) / dt;
    step = step==0 ? 1 : step;
    int count = 0;
    for (int h = begin_frame; h < end_frame; h+=step) {
        count++;
        if (log_flag == 1) {
            std::cout << " Time : " << h << std::endl;
            log_file << " ----------------------------------------------------- \n"
                     << " step " << std::setw(6) << h << std::endl
                     << " -------------------- " << std::endl;
        }

        std::vector<pinang::Vec3d> backbone1_nodes; // P atoms in the 1st backbone;
        std::vector<pinang::Vec3d> backbone1_dots;
        std::vector<pinang::Vec3d> backbone1_tangents; // Tangents at "P" atoms;
        std::vector<pinang::Vec3d> backbone1_normals; // Normal vectors at "P" atoms;

        std::vector<pinang::Vec3d> backbone2_nodes; // P atoms in the 2nd backbone;
        std::vector<pinang::Vec3d> backbone2_dots; // P atoms in the 2nd backbone;
        std::vector<pinang::Vec3d> backbone2_tangents;
        std::vector<pinang::Vec3d> backbone2_normals; // Normal vectors at "P" atoms;

        std::vector<pinang::Vec3d> base_positions1;
        std::vector<pinang::Vec3d> base_positions2;

        std::vector< std::vector<pinang::Vec3d> > genline1_lines; // Generating Lines!
        std::vector< std::vector<pinang::Vec3d> > genline2_lines;
        std::vector< std::vector<pinang::Vec3d> > genline1_line_tangents;
        std::vector< std::vector<pinang::Vec3d> > genline2_line_tangents;

        std::vector<pinang::Vec3d> axis_nodes;
        std::vector<pinang::Vec3d> axis_directions; // Axis direction vectors;

        std::vector<double> curv;

        pinang::Conformation confm = conformations[h];
        int len1 = phosphate_strand1_idx.size();
        int len2 = phosphate_strand2_idx.size();
        backbone1_nodes.clear();
        backbone2_nodes.clear();
        for (i = 0; i < len1; i++) {
            j = phosphate_strand1_idx[i];
            backbone1_nodes.push_back(confm.get_coor(j));
        }
        for (i = 0; i < len2; i++) {
            j = phosphate_strand2_idx[i];
            backbone2_nodes.push_back(confm.get_coor(j));
        }
        len1 = base_strand1_idx.size();
        len2 = base_strand2_idx.size();
        for (i = 0; i < len1; i++) {
            j = base_strand1_idx[i];
            base_positions1.push_back(confm.get_coor(j));
        }
        for (i = 0; i < len2; i++) {
            j = base_strand2_idx[i];
            base_positions2.push_back(confm.get_coor(j));
        }

        // -------------------- backbone --------------------
        gen_spline_fit(backbone1_nodes, 10, backbone1_dots, backbone1_tangents);
        gen_spline_fit(backbone2_nodes, 10, backbone2_dots, backbone2_tangents);

        // -------------------- genline --------------------
        genline1_lines.clear();
        genline2_lines.clear();
        genline1_line_tangents.clear();
        genline2_line_tangents.clear();
        for (i = 0; i < 10; i++) {
            std::vector<pinang::Vec3d> genline1_dots;  // Interpolation points;
            std::vector<pinang::Vec3d> genline1_nodes; // P atoms in the 1st backbone;
            std::vector<pinang::Vec3d> genline1_tangents; // Tangents at "P" atoms;
            std::vector<pinang::Vec3d> genline2_dots;
            std::vector<pinang::Vec3d> genline2_nodes; // P atoms in the 2nd backbone;
            std::vector<pinang::Vec3d> genline2_tangents;

            for (j = i; j < int(backbone1_nodes.size()); j += 10)
                genline1_nodes.push_back(backbone1_nodes[j]);
            for (j = i; j < int(backbone2_nodes.size()); j += 10)
                genline2_nodes.push_back(backbone2_nodes[j]);

            gen_spline_fit(genline1_nodes, 50, genline1_dots, genline1_tangents);
            gen_spline_fit(genline2_nodes, 50, genline2_dots, genline2_tangents);

            genline1_lines.push_back(genline1_dots);
            genline2_lines.push_back(genline2_dots);
            genline1_line_tangents.push_back(genline1_tangents);
            genline2_line_tangents.push_back(genline2_tangents);

            genline1_nodes.clear();
            genline2_nodes.clear();
            genline1_dots.clear();
            genline2_dots.clear();
            genline1_tangents.clear();
            genline2_tangents.clear();
        }

        // -------------------- backbone normal --------------------
        for (i = 0; i < int(backbone1_nodes.size()); i++) {
            int m = i / 10;
            int n = i % 10;
            pinang::Vec3d t1 = backbone1_tangents[i];
            pinang::Vec3d t2 = genline1_line_tangents[n][m];
            pinang::Vec3d nm = t2 ^ t1;
            pinang::Vec3d n0 = nm * (1.0 / nm.norm());
            backbone1_normals.push_back(n0);
        }
        for (i = 0; i < int(backbone2_nodes.size()); i++) {
            int m = i / 10;
            int n = i % 10;
            pinang::Vec3d t1 = backbone2_tangents[i];
            pinang::Vec3d t2 = genline2_line_tangents[n][m];
            pinang::Vec3d nm = t2 ^ t1;
            pinang::Vec3d n0 = nm * (1.0 / nm.norm());
            backbone2_normals.push_back(n0);
        }

        // -------------------- axis nodes --------------------
        for (i = 0; i < int(backbone1_normals.size()); i++) {
            pinang::Vec3d t1 = backbone1_normals[i];
            pinang::Vec3d t2 = genline1_line_tangents[i%10][i/10];
            pinang::Vec3d t3 = t1 ^ t2;
            pinang::Vec3d nm;       // normal vector of the fake plane
            double theta = 0;
            double alpha = 0;

            int cog_base_id = 0;
            double d_m = 1000.0;
            double _d_ = 0.0;
            for (j = 0; j < int(base_positions2.size()); j++) {
                _d_ = vec_distance(base_positions1[i], base_positions2[j]);
                if (_d_ < d_m) {
                    d_m = _d_;
                    cog_base_id = j;
                }
            }

            double circle_radius_min = 1000.0;
            pinang::Vec3d circle_center;
            for (theta = - angle_lim ; theta <= angle_lim; theta += pi_over_36) {
                nm = t2 + t3 * tan(theta);
                for (alpha = - angle_lim2 ; alpha <= angle_lim2; alpha += pi_over_60) {
                    nm = nm + t1 * tan(alpha);
                    nm = nm * (1.0 / nm.norm());
                    double _d0 = nm * backbone1_nodes[i]; // param in plane function!

                    double _dm = 0.0;   // min distances to the plane!
                    double _d  = 0.0;   // distances to the plane!
                    std::vector<pinang::Vec3d> intersects;
                    pinang::Vec3d insect;

                    for (j = 0; j < 10; j++) {
                        _dm = 100000000.0;
                        int s = (i - j) * 5;
                        int k_min = std::max(0, s - 50);
                        int k_max = std::min(s + 50, int(genline1_lines[j].size()));
                        for (int k = k_min; k < k_max; k++) {
                            _d = genline1_lines[j][k] * nm - _d0;
                            double d_tmp = _d > 0 ? _d : - _d;
                            if (d_tmp <= _dm) {
                                _dm = d_tmp;
                                insect = genline1_lines[j][k] - (nm * _d);
                            }
                        }
                        intersects.push_back(insect);
                    }
                    for (j = 0; j < 10; j++) {
                        _dm = 100000000.0;
                        int s = (cog_base_id - j) * 5;
                        int k_min = std::max(0, s - 50);
                        int k_max = std::min(s + 50, int(genline1_lines[j].size()));
                        for (int k = k_min; k < k_max; k++) {
                            _d = genline2_lines[j][k] * nm - _d0;
                            double d_tmp = _d > 0 ? _d : - _d;
                            if (d_tmp <= _dm) {
                                _dm = d_tmp;
                                insect = genline2_lines[j][k] - (nm * _d);
                            }
                        }
                        intersects.push_back(insect);
                    }

                    pinang::Vec3d com(0,0,0);
                    for (j = 0; j < int(intersects.size()); j++) {
                        com = com + intersects[j];
                    }
                    com = com * (1.0 / int(intersects.size()) );
                    double d_max = 0;   // largest distance from intersects to com
                    double d_tmp = 0;   // distance from intersects to com
                    for (j = 0; j < int(intersects.size()); j++) {
                        d_tmp = vec_distance(com, intersects[j]);
                        if (d_tmp > d_max)
                            d_max = d_tmp;
                    }
                    if (d_max < circle_radius_min) {
                        circle_radius_min = d_max;
                        circle_center = com;
                    }
                }
            }
            axis_nodes.push_back(circle_center);
        }

        // -------------------- axis direction --------------------
        for (i = 3; i < int(axis_nodes.size()) - 3; i++) {
            pinang::Vec3d s_tmp(0,0,0);
            pinang::Vec3d t_tmp(0,0,0);
            pinang::Vec3d v1(0,0,0);
            pinang::Vec3d v2(0,0,0);
            double t1 = 0;
            double t2 = 0;
            v1 = axis_nodes[i + 3] - axis_nodes[i];
            v2 = axis_nodes[i] - axis_nodes[i - 3];
            t1 = v1.sqr_norm();
            t2 = v2.sqr_norm();
            s_tmp = v1 * t2 + v2 * t1;
            t_tmp = s_tmp * (1 / s_tmp.norm());
            axis_directions.push_back(t_tmp);
        }

        // -------------------- curvature --------------------
        for (i = 4; i < int(axis_nodes.size())-4; i++) {
            j = i - 3;
            pinang::Vec3d v1 = axis_nodes[i + 1] - axis_nodes[i - 1];
            pinang::Vec3d v2 = axis_nodes[i + 2] - axis_nodes[i - 2];
            pinang::Vec3d v3 = axis_nodes[i + 3] - axis_nodes[i - 3];
            pinang::Vec3d dir1(0,0,0);
            pinang::Vec3d dir2(0,0,0);
            pinang::Vec3d dir3(0,0,0);
            double dist1 = v1.norm();
            double dist2 = v2.norm();
            double dist3 = v3.norm();
            double k1=0, k2=0, k3=0;
            double k_sum = 0;
            double k_ave = 0;
            int k_count = 0;
            int f_only_k1 = 0;
            if ( j >=1 && j < int(axis_directions.size()) - 1) {
                dir1 = axis_directions[j-1];
                dir2 = axis_directions[j+1];
                k1 = pinang::vec_angle(dir1, dir2) / dist1;
                if (dist1 > 3) {
                    k_sum += k1;
                    k_count += 1;
                    f_only_k1 = 1;
                }
            }
            if ( j >=2 && j < int(axis_directions.size()) - 2) {
                dir1 = axis_directions[j-2];
                dir2 = axis_directions[j+2];
                k2 = pinang::vec_angle(dir1, dir2) / dist2;
                if (dist2 > 3) {
                    k_sum += k2;
                    k_count += 1;
                    f_only_k1 = 0;
                }
            }
            if ( j >=3 && j < int(axis_directions.size()) - 3) {
                dir1 = axis_directions[j-3];
                dir2 = axis_directions[j+3];
                k3 = pinang::vec_angle(dir1, dir2) / dist3;
                if (dist3 > 3) {
                    k_sum += k3;
                    k_count += 1;
                    f_only_k1 = 0;
                }
            }
            if (k_count != 0 && f_only_k1 == 0) {
                k_ave = k_sum / k_count;
            }
            if (log_flag == 1) {
                log_file << std::setw(6) << i + 2 << "  "
                         << std::setw(10) << k_ave << "   "
                         << std::endl;
            }
            curv.push_back(k_ave);
        }
        curvature_all.push_back(curv);
        curv.clear();
    }

    std::cout << " Totally " << count << " frames analyzed."<< std::endl;

    std::vector<double> curv_ana;
    int len_strand = curvature_all[0].size();
    for (i = 0; i < len_strand; i++)
        curv_ana.push_back(0);
    for (i = 0; i < count; i++) {
        for (j = 0; j < len_strand; j++) {
            curv_ana[j] += curvature_all[i][j];
        }
    }
    for (i = 0; i < len_strand; i++) {
        curv_ana[i] /= count;
        crv_file << curv_ana[i] << std::endl;
    }



    // ending ------------------------------
    dcd_file.close();
    crv_file.close();
    log_file.close();

    return 0;
}



void gen_spline_fit(const std::vector<pinang::Vec3d>& nodes, int N,
                    std::vector<pinang::Vec3d>& dots,
                    std::vector<pinang::Vec3d>& tangents)
{
    int len1 = int(nodes.size());
    dots.clear();
    tangents.clear();

    pinang::Vec3d s_tmp(0,0,0);
    pinang::Vec3d t_tmp(0,0,0);
    pinang::Vec3d p_i(0,0,0);
    pinang::Vec3d b_r(0,0,0);
    pinang::Vec3d v1(0,0,0);
    pinang::Vec3d v2(0,0,0);
    double t1 = 0;
    double t2 = 0;
    double d = 0;
    double r = 0;

    for (int i = 0; i < len1; i++) {
        if (i == 0)
        {
            s_tmp = nodes[i + 1] - nodes[i];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            tangents.push_back(t_tmp);
        } else if (i == len1 - 1) {
            s_tmp = nodes[i] - nodes[i - 1];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            tangents.push_back(t_tmp);
        } else {
            v1 = nodes[i + 1] - nodes[i];
            v2 = nodes[i] - nodes[i - 1];
            t1 = v1.sqr_norm();
            t2 = v2.sqr_norm();
            s_tmp = v1 * t2 + v2 * t1;
            t_tmp = s_tmp * (1 / s_tmp.norm());
            tangents.push_back(t_tmp);
        }
    }

    // ~~~~~~~~~~ correction of the first and the last vectors ~~~~~~~~~~
    t1 = (tangents[1] * tangents[2]) * 2;
    v1 = tangents[1] * t1;
    tangents[0] = v1 - tangents[2]; // fixed tangent 0;

    t1 = (tangents[len1-3] * tangents[len1-2]) * 2;
    v1 = tangents[len1-2] * t1;
    tangents[len1-1] = v1 - tangents[len1-3]; // fixed last tangent;

    // ~~~~~~~~~~ spline fitting ~~~~~~~~~~
    for (int i = 0; i < len1 - 1; i++) {
        p_i = nodes[i];
        t_tmp = nodes[i+1] - nodes[i];
        d = t_tmp.norm();
        v1 = tangents[i];
        v2 = tangents[i+1];
        for (int j = 0; j < N; j++) {
            r = j * (1.0 / N);
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r ) + v2 * d * (r*r*r - r*r);
            dots.push_back(b_r);
        }
    }
    p_i = nodes[len1-1];
    dots.push_back(p_i);
}

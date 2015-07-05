#include "PDB.h"
#include "constants.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace std;

void gen_spline_fit(const std::vector<pinang::Vec3d>&, int N,
                    std::vector<pinang::Vec3d>&, std::vector<pinang::Vec3d>&);

int main(int argc, char *argv[])
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // whatever

    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~                  PINANG DNA curvature                  ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt, mod_index = 0;
    int mod_flag = 0;
    int in_flag = 0;
    int inp_flag = 0;

    std::string infilename = "some.pdb";

    std::string inp_name = "curve.inp";
    std::string axis_name = "_axis.pdb";
    std::string norm_name = "_norm.pdb";
    std::string groove_name = "_groove.pdb";
    std::string back_name = "_backbone.pdb";
    std::string gline_name = "_generating_lines.pdb";
    std::string out_name = "_curve.dat";

    while ((opt = getopt(argc, argv, "o:x:b:g:i:m:f:h")) != -1) {
        switch (opt) {
        case 'o':
            out_name = optarg;
            break;
        case 'i':
            inp_name = optarg;
            inp_flag = 1;
            break;
        case 'x':
            axis_name = optarg;
            break;
        case 'b':
            back_name = optarg;
            break;
        case 'g':
            gline_name = optarg;
            break;
        case 'm':
            mod_index = atoi(optarg);
            mod_flag = 1;
            break;
        case 'f':
            infilename = optarg;
            in_flag = 1;
            break;
        case 'h':
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                      << " [-b _backbone.pdb] [-g _generating_lines.pdb]"
                      << " [-m module] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                      << " [-b _backbone.pdb] [-g _generating_lines.pdb]"
                      << " [-m module] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag || !inp_flag)
    {
        std::cout << " ERROR: need parameter for option -f and -i: "
                  << std::endl << " Usage: " << argv[0]
                  << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                  << " [-b _backbone.pdb] [-g _generating_lines.pdb]"
                  << " [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    std::ifstream inp_file(inp_name.c_str());
    std::ofstream back_file(back_name.c_str());
    std::ofstream gline_file(gline_name.c_str());
    std::ofstream axis_file(axis_name.c_str());
    std::ofstream norm_file(norm_name.c_str());
    std::ofstream groove_file(groove_name.c_str());
    std::ofstream out_file(out_name.c_str());

    if (mod_flag != 1) {
        if (pdb1.n_models() == 1)
        {
            mod_index = 1;
        } else {
            std::cout << " Please choose a MODULE: " ;
            std::cin >> mod_index;
        }
    }

    std::cout << " Analyzing DNA curvature of MODULE " << mod_index
              << " of " << infilename  << " ... " << std::endl
              << std::endl;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // main
    //  _ _                 _      _ _
    // ( | )_ __ ___   __ _(_)_ __( | )
    //  V V| '_ ` _ \ / _` | | '_ \V V
    //     | | | | | | (_| | | | | |
    //     |_| |_| |_|\__,_|_|_| |_|
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    std::vector<pinang::Vec3d> backbone1_dots;  // Interpolation points;
    std::vector<pinang::Vec3d> backbone1_nodes; // P atoms in the 1st backbone;
    std::vector<pinang::Vec3d> backbone1_tangents; // Tangents at "P" atoms;
    std::vector<pinang::Vec3d> backbone1_normals; // Normal vectors at "P" atoms;

    std::vector<pinang::Vec3d> backbone2_dots;
    std::vector<pinang::Vec3d> backbone2_nodes; // P atoms in the 2nd backbone;
    std::vector<pinang::Vec3d> backbone2_tangents;
    std::vector<pinang::Vec3d> backbone2_normals; // Normal vectors at "P" atoms;

    std::vector< std::vector<pinang::Vec3d> > genline1_lines; // Generating Lines!
    std::vector< std::vector<pinang::Vec3d> > genline2_lines;
    std::vector< std::vector<pinang::Vec3d> > genline1_line_tangents;
    std::vector< std::vector<pinang::Vec3d> > genline2_line_tangents;

    std::vector<pinang::Vec3d> base_positions1;
    std::vector<pinang::Vec3d> base_positions2;

    std::vector<pinang::Vec3d> axis_dots;
    std::vector<pinang::Vec3d> axis_nodes;
    std::vector<pinang::Vec3d> axis_directions; // Axis direction vectors;
    std::vector<double> helix_width;

    std::vector<double> base_rise;
    std::vector<double> bases_per_turn;

    std::vector<double> major_groove_width;
    std::vector<double> minor_groove_width;



    /* ============================================================
    //      _                   _
    //     (_)_ __  _ __  _   _| |_
    //     | | '_ \| '_ \| | | | __|
    //     | | | | | |_) | |_| | |_
    //     |_|_| |_| .__/ \__,_|\__|
    //              |_|
    // ============================================================
    */
    std::cout << " 1. Read in structures ..." << std::endl;
    int i = 0;
    int j = 0;

    int flg_grp_1 = 0;
    int flg_grp_2 = 0;
    std::vector<char> chain1_id;
    std::vector<char> chain2_id;
    std::string inp_line;
    std::string tmp_str;
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
            for (i = 0; i < int(strs.size()); i++) {
                tmp_s = strs[i];
                char tmp_c = 0;
                // tmp_sstr.str(tmp_s);
                tmp_sstr.str(strs[i]);
                tmp_sstr >> tmp_c;
                chain1_id.push_back(tmp_c);
                tmp_sstr.clear();
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
            for (i = 0; i < int(strs.size()); i++) {
                tmp_s = strs[i];
                char tmp_c = 0;
                // tmp_sstr.str(tmp_s);
                tmp_sstr.str(strs[i]);
                tmp_sstr >> tmp_c;
                if (tmp_c == ' ' or tmp_c == 0)
                    continue;
                chain2_id.push_back(tmp_c);
                tmp_sstr.clear();
            }
        }
    }
    if (flg_grp_1 == 0 || chain1_id.size() == 0)
    {
        std::cout << " ERROR: STRAND1 not found!" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (flg_grp_2 == 0 || chain2_id.size() == 0)
    {
        std::cout << " ERROR: STRAND2 not found!" << std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout << "Strand 1: ";
    // for (i = 0; i < int(chain1_id.size()); i++) {
    //     std::cout << chain1_id[i] << ", ";
    // }
    // std::cout << std::endl;
    // std::cout << "Strand 2: ";
    // for (i = 0; i < int(chain2_id.size()); i++) {
    //     std::cout << chain2_id[i] << ", ";
    // }
    // std::cout << std::endl;
    inp_file.close();

    pinang::Model mdl0 = pdb1.m_model(mod_index-1);
    pinang::Chain chain_tmp;
    int mdl_size = mdl0.m_model_size();
    for (i = 0; i < int(chain1_id.size()); i++)
        for (j = 0; j < mdl_size; ++j)
            if (mdl0.m_chain(j).chain_ID() == chain1_id[i]) {
                chain_tmp = mdl0.m_chain(j);
                int len1 = chain_tmp.m_chain_length();
                for (int k = 1; k < len1; k++) { // start from 1! because residue 0 has no P!
                    backbone1_nodes.push_back(chain_tmp.m_residue(k).m_P().coordinates());
                }
                for (int k = 0; k < len1; k++) {
                    base_positions1.push_back(chain_tmp.m_residue(k).m_B().coordinates());
                }
                chain_tmp.reset();
            }
    for (i = 0; i < int(chain2_id.size()); i++)
        for (j = 0; j < mdl_size; ++j)
            if (mdl0.m_chain(j).chain_ID() == chain2_id[i]) {
                chain_tmp = mdl0.m_chain(j);
                int len1 = chain_tmp.m_chain_length();
                for (int k = 1; k < len1; k++) { // start from 1! because residue 0 has no P!
                    backbone2_nodes.push_back(chain_tmp.m_residue(k).m_P().coordinates());
                }
                for (int k = 0; k < len1; k++) {
                    base_positions2.push_back(chain_tmp.m_residue(k).m_B().coordinates());
                }
                chain_tmp.reset();
            }
    if (backbone1_nodes.size() == 0 || base_positions1.size() == 0)
    {
        std::cout << " ERROR: STRAND1 read in error!" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (backbone2_nodes.size() == 0 || base_positions2.size() == 0)
    {
        std::cout << " ERROR: STRAND2 read in error!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << " ... done." << std::endl;


    /* ============================================================
    //  _                _    _
    // | |__   __ _  ___| | _| |__   ___  _ __   ___
    // | '_ \ / _` |/ __| |/ / '_ \ / _ \| '_ \ / _ \
    // | |_) | (_| | (__|   <| |_) | (_) | | | |  __/
    // |_.__/ \__,_|\___|_|\_\_.__/ \___/|_| |_|\___|
    // ============================================================
    */
    std::cout << " 2. Generating backbone curves ..." << std::endl;
    gen_spline_fit(backbone1_nodes, 10, backbone1_dots, backbone1_tangents);
    gen_spline_fit(backbone2_nodes, 10, backbone2_dots, backbone2_tangents);

    // ============================ Output to PDB ============================
    int k = backbone1_dots.size();
    int l = backbone1_nodes.size();
    for (i = 0; i < int(backbone1_dots.size()); i++) {
        if (i % 10 == 0) {
            back_file << std::setw(6) << "HETATM" << std::setw(5) << i+1 << " "
                      << std::setw(4) << "C   "   << std::setw(1) << " "
                      << std::setw(3) << "CUR" << " " << std::setw(1) << "A"
                      << std::setw(4) << i/10+2 << std::setw(1) << " " << "   "
                      << backbone1_dots[i] << std::endl;
        } else {
            back_file << std::setw(6) << "HETATM" << std::setw(5) << i+1 << " "
                      << std::setw(4) << "O   " << std::setw(1) << " "
                      << std::setw(3) << "CUR" << " " << std::setw(1) << "A"
                      << std::setw(4) << i/10+2 << std::setw(1) << " " << "   "
                      << backbone1_dots[i] << std::endl;
        }
    }
    back_file << "TER" << std::endl;
    for (i = 0; i < int(backbone2_dots.size()); i++) {
        if (i % 10 == 0) {
            back_file << std::setw(6) << "HETATM" << std::setw(5) << i+1+k << " "
                      << std::setw(4) << "C   "   << std::setw(1) << " "
                      << std::setw(3) << "CUR" << " " << std::setw(1) << "B"
                      << std::setw(4) << i/10+l+3 << std::setw(1) << " " << "   "
                      << backbone2_dots[i] << std::endl;
        } else {
            back_file << std::setw(6) << "HETATM" << std::setw(5) << i+1+k << " "
                      << std::setw(4) << "O   " << std::setw(1) << " "
                      << std::setw(3) << "CUR" << " " << std::setw(1) << "B"
                      << std::setw(4) << i/10+l+3 << std::setw(1) << " " << "   "
                      << backbone2_dots[i] << std::endl;
        }
    }

    // connecting points!
    for (i = 0; i < int(backbone1_dots.size())-1; i++) {
        back_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1
                   << std::setw(5) << i+2
                   << std::endl;
    }
    for (i = 0; i < int(backbone2_dots.size())-1; i++) {
        back_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1+k
                   << std::setw(5) << i+2+k
                   << std::endl;
    }

    std::cout << " ... done." << std::endl;

    /* ============================================================
    //                     _ _
    //   __ _  ___ _ __   | (_)_ __   ___
    //  / _` |/ _ \ '_ \  | | | '_ \ / _ \
    // | (_| |  __/ | | | | | | | | |  __/
    //  \__, |\___|_| |_| |_|_|_| |_|\___|
    //  |___/
    // ============================================================
    */
    std::cout << " 3. Generating lines of cylinders ..." << std::endl;
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

        // ============================ Output to PDB ==========================
        int k = int(genline1_dots.size());
        int l = int(backbone1_nodes.size());
        for (j = 0; j < k; j++) {
            if (j % 50 == 0) {
                gline_file << std::setw(6) << "HETATM" << std::setw(5) << j+1 << " "
                           << std::setw(4) << "C   "   << std::setw(1) << " "
                           << std::setw(3) << "GEN" << " " << std::setw(1) << "A"
                           << std::setw(4) << (j/50) * 10 + i + 1 << std::setw(1) << " " << "   "
                           << genline1_dots[j] << std::endl;
            } else {
                gline_file << std::setw(6) << "HETATM" << std::setw(5) << j+1 << " "
                           << std::setw(4) << "O   " << std::setw(1) << " "
                           << std::setw(3) << "GEN" << " " << std::setw(1) << "A"
                           << std::setw(4) << (j/50) * 10 + i + 1 << std::setw(1) << " " << "   "
                           << genline1_dots[j] << std::endl;
            }
        }
        gline_file << "TER" << std::endl;
        for (j = 0; j < int(genline2_dots.size()); j++) {
            if (j % 10 == 0) {
                gline_file << std::setw(6) << "HETATM" << std::setw(5) << j+1+k << " "
                           << std::setw(4) << "C   "   << std::setw(1) << " "
                           << std::setw(3) << "GEN" << " " << std::setw(1) << "B"
                           << std::setw(4) << j/50 * 10 + l + i + 1 << std::setw(1) << " " << "   "
                           << genline2_dots[j] << std::endl;
            } else {
                gline_file << std::setw(6) << "HETATM" << std::setw(5) << j+1+k << " "
                           << std::setw(4) << "O   " << std::setw(1) << " "
                           << std::setw(3) << "GEN" << " " << std::setw(1) << "B"
                           << std::setw(4) << j/50 * 10 + l + i + 1 << std::setw(1) << " " << "   "
                           << genline2_dots[j] << std::endl;
            }
        }
        gline_file << "TER" << std::endl;

        // connecting points!
        for (j = 0; j < int(genline1_dots.size())-1; j++) {
            gline_file << std::setw(6) << "CONECT"
                      << std::setw(5) << j + 1
                      << std::setw(5) << j + 2
                      << std::endl;
        }
        for (j = 0; j < int(genline2_dots.size())-1; j++) {
            gline_file << std::setw(6) << "CONECT"
                      << std::setw(5) << j+1+k
                      << std::setw(5) << j+2+k
                      << std::endl;
        }
        genline1_nodes.clear();
        genline2_nodes.clear();
        genline1_dots.clear();
        genline2_dots.clear();
        genline1_tangents.clear();
        genline2_tangents.clear();
    }
    std::cout << " ... done." << std::endl;

    /* ==================================================================
    //  _          _ _            _ _               _   _
    // | |__   ___| (_)_  __   __| (_)_ __ ___  ___| |_(_) ___  _ __
    // | '_ \ / _ \ | \ \/ /  / _` | | '__/ _ \/ __| __| |/ _ \| '_ \
    // | | | |  __/ | |>  <  | (_| | | | |  __/ (__| |_| | (_) | | | |
    // |_| |_|\___|_|_/_/\_\  \__,_|_|_|  \___|\___|\__|_|\___/|_| |_|
    // ==================================================================
    */
    std::cout << " 4. Calculating helix axis ..." << std::endl;
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
    // -------------------- output normal vectors for P --------------------
    for (j = 0; j < int(backbone1_nodes.size()); j++) {
        norm_file << std::setw(6) << "HETATM" << std::setw(5) << 3 * j+1 << " "
                  << std::setw(4) << "C   " << std::setw(1) << " "
                  << std::setw(3) << "AXS" << " " << std::setw(1) << "A"
                  << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                  << backbone1_nodes[j] << std::endl;
        norm_file << std::setw(6) << "HETATM" << std::setw(5) << 3 * j+2 << " "
                  << std::setw(4) << "O   " << std::setw(1) << " "
                  << std::setw(3) << "AXS" << " " << std::setw(1) << "A"
                  << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                  << backbone1_nodes[j] + backbone1_normals[j] * 5 << std::endl;
    }
    // connecting points!
    for (j = 0; j < int(backbone1_nodes.size()); j++) {
        norm_file << std::setw(6) << "CONECT"
                   << std::setw(5) << 3 * j + 1
                   << std::setw(5) << 3 * j + 2
                   << std::endl;
    }


    // ============= find the helix center and tangent for every P =============
    double angle_lim = pinang::g_pi / 3;   // 60 degree;
    double angle_lim2 = pinang::g_pi / 15;   // 12 degree;
    double pi_over_36 = pinang::g_pi / 36; // 5 degree;
    double pi_over_60 = pinang::g_pi / 60; // 3 degree;
    for (i = 0; i < int(backbone1_normals.size()); i++) {
        // if (i != 52 && i != 53 && i != 51)
        // if (i != 52)
        //     continue;

        pinang::Vec3d t1 = backbone1_normals[i];
        pinang::Vec3d t2 = genline1_line_tangents[i%10][i/10];
        pinang::Vec3d t3 = t1 ^ t2;
        // t3 = t3 * (1.0 / t3.norm());
        pinang::Vec3d nm;       // normal vector of the fake plane
        double theta = 0;
        double alpha = 0;

        // find closest base on genline2...
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
        // std::cout << i << "  vs  " << cog_base_id << " -> " << i + cog_base_id << std::endl;

        double circle_radius_min = 1000.0;
        pinang::Vec3d circle_center;
        for (theta = - angle_lim ; theta <= angle_lim; theta += pi_over_36) {
            nm = t2 + t3 * tan(theta);
            for (alpha = - angle_lim2 ; alpha <= angle_lim2; alpha += pi_over_60) {
                nm = nm + t1 * tan(alpha);
                nm = nm * (1.0 / nm.norm());
                // STEP 1: calculate the plane! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // Plane: a x + b y + c z = _d0
                double _d0 = nm * backbone1_nodes[i]; // param in plane function!

                // STEP 2: get intersections for each generating line ~~~~~~~~~~~~~~
                double _dm = 0.0;   // min distances to the plane!
                double _d  = 0.0;   // distances to the plane!
                std::vector<pinang::Vec3d> intersects;
                pinang::Vec3d insect;

                for (j = 0; j < int(genline1_lines.size()); j++) {
                    _dm = 100000000.0;
                    // find the closest points on genlines...
                    int s = (i - j) * 5;
                    int k_min = std::max(0, s - 50);
                    int k_max = std::min(s + 50, int(genline1_lines[j].size()));
                    // find the closets point and the crspnding intersect
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
                for (j = 0; j < int(genline2_lines.size()); j++) {
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
                // STEP 3: find out the minimal circle cover all dots ~~~~~~~~~~~~~~
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
                // STEP 4: return the helix center and helix width!  ~~~~~~~~~~~~~~~
                if (d_max < circle_radius_min) {
                    circle_radius_min = d_max;
                    circle_center = com;
                }
            }
        }
        axis_nodes.push_back(circle_center);
        helix_width.push_back(circle_radius_min);
    }

    // Calculate Helix DIRECTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    // Output to PDB...
    for (j = 4; j < int(axis_nodes.size())-4; j++) {
        axis_file << std::setw(6) << "HETATM" << std::setw(5) << j+1 << " "
                  << std::setw(4) << "C   " << std::setw(1) << " "
                  << std::setw(3) << "AXS" << " " << std::setw(1) << "A"
                  << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                  << axis_nodes[j] << std::endl;
        norm_file << std::setw(6) << "HETATM" << std::setw(5) << 3 * j+3 << " "
                  << std::setw(4) << "C   " << std::setw(1) << " "
                  << std::setw(3) << "AXS" << " " << std::setw(1) << "A"
                  << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                  << axis_nodes[j] << std::endl;
    }
    // connecting points!
    for (j = 4; j < int(axis_nodes.size())-4; j++) {
        axis_file << std::setw(6) << "CONECT"
                   << std::setw(5) << j + 1
                   << std::setw(5) << j + 2
                   << std::endl;
        norm_file << std::setw(6) << "CONECT"
                  << std::setw(5) << 3* j + 2
                  << std::setw(5) << 3* j + 3
                  << std::endl;
    }

    // ============================================================
    // Calculate Curvature!
    out_file << "------------------------------------------------------------"
             << std::endl;
    out_file << "DNA curvature: " << std::endl
             << std::setw(6) << "id" << "  "
             << std::setw(10) << "k1" << "   "
             << std::setw(10) << "k2" << "   "
             << std::setw(10) << "k3" << "   "
             << std::setw(10) << "k_ave" << "   "
             << std::endl;
    out_file << "------------------------------------------------------------"
             << std::endl;

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
        // ----- i - 1 : i + 1
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

        out_file << std::setw(6) << i + 2 << "  "
                 << std::setw(10) << k1 << "   "
                 << std::setw(10) << k2 << "   "
                 << std::setw(10) << k3 << "   "
                 << std::setw(10) << k_ave << "   "
                 << std::endl;
    }
    out_file << "============================================================\n\n\n"
             << std::endl;

    // ---------- Base rise ----------
    // pinang::Vec3d base_delta;
    // double proj = 0;
    // for (j = 4; j < int(base_positions1.size())-5; j++) {
    //     base_delta = base_positions1[j+1] - base_positions1[j];
    //     proj = base_delta * axis_directions[j-4];
    //     base_rise.push_back(abs(proj));
    //     // std::cout << j << " - " << j + 1 << " : "
    //     //           << base_delta.norm() << "    "
    //     //           << abs(proj) << std::endl;
    // }
    std::cout << " ... done." << std::endl;

    /* =========================================================================
    //   ____                                     _     _ _   _
    //  / ___|_ __ ___   _____   _____  __      _(_) __| | |_| |__
    // | |  _| '__/ _ \ / _ \ \ / / _ \ \ \ /\ / / |/ _` | __| '_ \
    // | |_| | | | (_) | (_) \ V /  __/  \ V  V /| | (_| | |_| | | |
    //  \____|_|  \___/ \___/ \_/ \___|   \_/\_/ |_|\__,_|\__|_| |_|
    // =========================================================================
    */
    std::cout << " 5. Calculating groove width ..." << std::endl;
    out_file << "------------------------------------------------------------"
             << std::endl;
    out_file << std::setw(6) << "id" << "   "
            << std::setw(20) << "minor groove width" << "   "
            << std::setw(20) <<  "major groove width" << std::endl;
    out_file << "------------------------------------------------------------"
             << std::endl;

    for (i = 0; i < int(axis_directions.size())-1; i++) {
        // step 1: plane perpendicular to direction Ox, O is the current point -
        // Plane: a x + b y + c z = _d0
        j = i + 4;
        pinang::Vec3d nm = axis_directions[i];
        pinang::Vec3d _O = axis_nodes[j];
        double _d0 = nm * _O; // param in plane function!
        double _d = 0;
        // step 2: intersections with backbones it1 and it2 --------------------
        pinang::Vec3d it1;
        pinang::Vec3d it2;
        // backbone1
        double _dm = 100000000.0;
        int k_min = std::max(0, j * 10 - 50);
        int k_max = std::min(j * 10 + 50, int(backbone1_dots.size()));
        int kb1 = 0;
        int kb2 = 0;
        for (int k = k_min; k < k_max; k++) {
            _d = backbone1_dots[k] * nm - _d0;
            double d_tmp = _d > 0 ? _d : - _d;
            if (d_tmp <= _dm) {
                _dm = d_tmp;
                it1 = backbone1_dots[k] - (nm * _d);
                kb1 = k;
            }
        }
        // backbone2
        _dm = 100000000.0;
        k = int(backbone2_nodes.size());
        k_min = std::max(0, (k - j) * 10 - 100);
        k_max = std::min((k - j) * 10 + 100, int(backbone2_dots.size()));
        for (int k = k_min; k < k_max; k++) {
            _d = backbone2_dots[k] * nm - _d0;
            double d_tmp = _d > 0 ? _d : - _d;
            if (d_tmp <= _dm) {
                _dm = d_tmp;
                it2 = backbone2_dots[k] - (nm * _d);
                kb2 = k;
            }
        }

        // step 3: find a point in groove --------------------------------------
        pinang::Vec3d vtmp;
        pinang::Vec3d v1;
        pinang::Vec3d v2;
        pinang::Vec3d groove_D1;
        pinang::Vec3d nm_plane;

        pinang::Vec3d I1_0;   // backbone 1
        pinang::Vec3d I3_0;   // backbone 1

        pinang::Vec3d I2_0;   // backbone 2
        pinang::Vec3d I4_0;   // backbone 2

        vtmp = it1 - _O;
        v1 = vtmp * (1.0 / vtmp.norm());
        vtmp = it2 - _O;
        v2 = vtmp * (1.0 / vtmp.norm());
        vtmp = v1 + v2;
        groove_D1 = nm ^ v1;    // ****************************** OMG~
        if (vtmp.norm() > 0.05) {
            if (vtmp * groove_D1 > 0)
                groove_D1 = vtmp * (1.0 / vtmp.norm());
            else
                groove_D1 = vtmp * ( - 1.0 / vtmp.norm());
        }
        // ============================================================
        // !!! KEY~ !!!
        // step 4: rotate the test plane around groove_D
        double minor_g_w = 100000.0;   // local minor groove width----
        double major_g_w = 100000.0;
        pinang::Vec3d t3 = groove_D1 ^ nm; // plane normal vector
        double theta = 0;
        double angle_lim = pinang::g_pi / 2;   // 90 degree;
        for (theta = - angle_lim; theta <= 0; theta += pi_over_60) {
            pinang::Vec3d I1;   // backbone 1
            pinang::Vec3d I2;   // backbone 2

            // two intersection lines' I1-I2 and I3-I4

            nm_plane = t3 + nm * tan(theta);
            nm_plane = nm_plane * (1.0 / nm_plane.norm());
            _d0 = nm_plane * _O; // param in plane function!

            // ------------------------------ backbone 1 intersects
            _dm = 1000000.0;
            k_min = kb1;
            k_max = std::min(kb1 + 100, int(backbone1_dots.size()));
            for (k = k_min; k < k_max; k++) {
                _d = backbone1_dots[k] * nm_plane - _d0;
                double d_tmp = _d > 0 ? _d : - _d;
                if (d_tmp <= _dm) {
                    _dm = d_tmp;
                    I1 = backbone1_dots[k];
                }
            }

            // ------------------------------ backbone 2 intersects
            _dm = 1000000.0;
            k_min = kb2;
            k_max = std::min(kb2 + 100, int(backbone2_dots.size()));
            for (k = k_min; k < k_max; k++) {
                _d = backbone2_dots[k] * nm_plane - _d0;
                double d_tmp = _d > 0 ? _d : - _d;
                if (d_tmp <= _dm) {
                    _dm = d_tmp;
                    I2 = backbone2_dots[k];
                }
            }

            pinang::Vec3d min_g_v = I1 - I2;
            double mingw = min_g_v.norm();
            if (mingw < minor_g_w) {
                minor_g_w = mingw;
                I1_0 = I1;
                I2_0 = I2;
            }
        }
        double ang1 = acos(nm * backbone1_tangents[kb1/10]);
        double ang2 = acos(nm * backbone2_tangents[kb2/10]);
        if (ang1 > 1.570795)
            ang1 = pinang::g_pi - ang1;
        if (ang2 > 1.570795)
            ang2 = pinang::g_pi - ang2;
        angle_lim = pinang::g_pi / 2 - min(ang1, ang2);
        for (theta = 0; theta <= angle_lim; theta += pi_over_60) {
            pinang::Vec3d I3;   // backbone 1
            pinang::Vec3d I4;   // backbone 2

            nm_plane = t3 + nm * tan(theta);
            nm_plane = nm_plane * (1.0 / nm_plane.norm());
            _d0 = nm_plane * _O; // param in plane function!

            // ------------------------------ backbone 1 intersects
            _dm = 1000000.0;
            k_min = std::max(0, kb1 - 100);
            k_max = kb1;
            for (k = k_max; k > k_min; k--) {
                _d = backbone1_dots[k] * nm_plane - _d0;
                double d_tmp = _d > 0 ? _d : - _d;
                if (d_tmp <= _dm) {
                    _dm = d_tmp;
                    I3 = backbone1_dots[k];
                }
            }

            // ------------------------------ backbone 2 intersects
            _dm = 1000000.0;
            k_min = std::max(0, kb2 - 100);
            k_max = kb2;
            for (k = k_max; k > k_min; k--) {
                _d = backbone2_dots[k] * nm_plane - _d0;
                double d_tmp = _d > 0 ? _d : - _d;
                if (d_tmp <= _dm) {
                    _dm = d_tmp;
                    I4 = backbone2_dots[k];
                }
            }
            // if ( (I4 - _O) * groove_D1 > 0 )
            //     continue;

            pinang::Vec3d maj_g_v = I3 - I4;
            double majgw = maj_g_v.norm();
            if (majgw < major_g_w) {
                major_g_w = majgw;
                I3_0 = I3;
                I4_0 = I4;
            }
        }
        out_file << std::setw(6) << i+6 << "   "
                << std::setw(20) << minor_g_w << "   "
                << std::setw(20) <<  major_g_w << std::endl;
        if (minor_g_w < 25) {
            groove_file << std::setw(6) << "HETATM" << std::setw(5) << 4 * i+1 << " "
                        << std::setw(4) << "C   " << std::setw(1) << " "
                        << std::setw(3) << "GRV" << " " << std::setw(1) << "A"
                        << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                        << I1_0 << std::endl;
            groove_file << std::setw(6) << "HETATM" << std::setw(5) << 4 * i+2 << " "
                        << std::setw(4) << "N   " << std::setw(1) << " "
                        << std::setw(3) << "GRV" << " " << std::setw(1) << "B"
                        << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                        << I2_0 << std::endl;
            groove_file << std::setw(6) << "CONECT"
                        << std::setw(5) << 4*i + 1
                        << std::setw(5) << 4*i + 2
                        << std::endl;
            minor_groove_width.push_back(minor_g_w);
        }
        if (major_g_w < 25) {
            groove_file << std::setw(6) << "HETATM" << std::setw(5) << 4 * i+3 << " "
                        << std::setw(4) << "O   " << std::setw(1) << " "
                        << std::setw(3) << "GRV" << " " << std::setw(1) << "A"
                        << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                        << I3_0 << std::endl;
            groove_file << std::setw(6) << "HETATM" << std::setw(5) << 4 * i+4 << " "
                        << std::setw(4) << "S   " << std::setw(1) << " "
                        << std::setw(3) << "GRV" << " " << std::setw(1) << "B"
                        << std::setw(4) << j+1 << std::setw(1) << " " << "   "
                        << I4_0 << std::endl;
            groove_file << std::setw(6) << "CONECT"
                        << std::setw(5) << 4*i + 3
                        << std::setw(5) << 4*i + 4
                        << std::endl;
            major_groove_width.push_back(major_g_w);
        }
    }

    std::cout << " ... done." << std::endl;


    // ----------------------------------------------------------------------
    back_file.close();
    gline_file.close();
    axis_file.close();
    norm_file.close();
    groove_file.close();
    out_file.close();

    return 0;
}





/* FUNCTION gen_spline_fit
//                                 _ _                 __ _ _
//   __ _  ___ _ __      ___ _ __ | (_)_ __   ___     / _(_) |_
//  / _` |/ _ \ '_ \    / __| '_ \| | | '_ \ / _ \   | |_| | __|
// | (_| |  __/ | | |   \__ \ |_) | | | | | |  __/   |  _| | |_
//  \__, |\___|_| |_|___|___/ .__/|_|_|_| |_|\___|___|_| |_|\__|
//  |___/          |_____|  |_|                 |_____|
*/
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

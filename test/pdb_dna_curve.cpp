#include "PDB.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>

using namespace std;

void gen_spline_fit(const std::vector<pinang::Vec3d>&,
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

    std::string infilename = "some.pdb";

    std::string axis_name = "_axis.pdb";
    std::string back_name = "_backbone.pdb";
    std::string out_name = "_curve.dat";

    while ((opt = getopt(argc, argv, "o:x:b:m:f:h")) != -1) {
        switch (opt) {
        case 'o':
            out_name = optarg;
            break;
        case 'x':
            axis_name = optarg;
            break;
        case 'b':
            back_name = optarg;
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
                      << " [-b _backbone.pdb] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " Usage: "
                      << argv[0]
                      << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                      << " [-b _backbone.pdb] [-m module] [-h]"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!in_flag)
    {
        std::cout << " ERROR: need parameter for option -f: " << std::endl
                  << " Usage: "
                  << argv[0]
                  << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
                  << " [-b _backbone.pdb] [-m module] [-h]"
                  << std::endl;
        exit(EXIT_SUCCESS);
    }
    pinang::PDB pdb1(infilename);

    std::ofstream back_file(back_name.c_str());
    std::ofstream axis_file(axis_name.c_str());
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

    std::vector<pinang::Vec3d> backbone2_dots;
    std::vector<pinang::Vec3d> backbone2_nodes; // P atoms in the 2nd backbone;
    std::vector<pinang::Vec3d> backbone2_tangents;

    std::vector<pinang::Vec3d> base_positions;

    std::vector<pinang::Vec3d> axis_dots;
    std::vector<pinang::Vec3d> axis_nodes;
    std::vector<pinang::Vec3d> axis_directions; // Axis direction vectors;

    std::vector<double> base_rise;
    std::vector<double> helix_width;
    std::vector<double> bases_per_turn;

    std::vector<double> major_groove_width;
    std::vector<double> minor_groove_width;

    int i  = 0;
    pinang::Chain c1 = pdb1.m_model(mod_index-1).m_chain(0);
    pinang::Chain c2 = pdb1.m_model(mod_index-1).m_chain(1);
    int len1 = c1.m_chain_length();
    int len2 = c2.m_chain_length();
    for (i = 1; i < len1; i++) { // start from 1! because residue 0 has no P!
        backbone1_nodes.push_back(c1.m_residue(i).m_P().coordinates());
    }
    for (i = 1; i < len2; i++) {
        backbone2_nodes.push_back(c2.m_residue(i).m_P().coordinates());
    }
    for (i = 0; i < len1; i++) {
        base_positions.push_back(c1.m_residue(i).m_B().coordinates());
    }

    // ==================== backbone 1 calculation ====================
    len1 = int(backbone1_nodes.size());
    len2 = int(backbone2_nodes.size());

    gen_spline_fit(backbone1_nodes, backbone1_dots, backbone1_tangents);

    // ==================== backbone 2 calculation ====================
    gen_spline_fit(backbone2_nodes, backbone2_dots, backbone2_tangents);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ==================== axis calculation ====================


    // ============================ Output to PDB ============================
    // -------------------- backbone.pdb --------------------
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



    back_file.close();
    axis_file.close();
    out_file.close();

    return 0;
}

void gen_spline_fit(const std::vector<pinang::Vec3d>& nodes,
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
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r ) + v2 * d * (r*r*r - r*r);
            dots.push_back(b_r);
        }
    }
    p_i = nodes[len1-1];
    dots.push_back(p_i);
}

#include "PDB.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

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
    std::cout << " Usage: "
              << argv[0]
              << " -f some.pdb [-o _curve.dat] [-x _axis.pdb] \n"
              << " [-b _backbone.pdb] [-m module] [-h]"
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
              << " of " << infilename  << " ... "
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

    const double p_D_phosphate = 17.45;

    std::vector<pinang::Vec3d> curve1_dots;  // Interpolation points;
    std::vector<pinang::Vec3d> curve1_nodes; // "P" atoms in the 1st backbone;
    std::vector<pinang::Vec3d> curve1_tangents; // Tangents at "P" atoms;

    std::vector<pinang::Vec3d> curve2_dots;
    std::vector<pinang::Vec3d> curve2_nodes; // "P " atoms in the 2nd backbone;
    std::vector<pinang::Vec3d> curve2_tangents;

    std::vector<pinang::Vec3d> cross_vecs;
    std::vector<pinang::Vec3d> base_positions;

    std::vector<pinang::Vec3d> axis_dots;
    std::vector<pinang::Vec3d> axis_nodes;
    std::vector<pinang::Vec3d> axis_directions; // Axis direction vectors;

    std::vector<double> base_rise;
    std::vector<double> helix_width;
    std::vector<double> bases_per_turn;

    std::vector<double> major_groove_width;
    std::vector<double> minor_groove_width;

    double d = 0;               // tmp for base_rise;
    double r = 0;               // tmp for helix width;

    int i  = 0;
    pinang::Chain c1 = pdb1.m_model(mod_index-1).m_chain(0);
    pinang::Chain c2 = pdb1.m_model(mod_index-1).m_chain(1);
    int len1 = c1.m_chain_length();
    int len2 = c2.m_chain_length();
    for (i = 1; i < len1; i++) {
        curve1_nodes.push_back(c1.m_residue(i).m_P().coordinates());
    }
    for (i = 1; i < len2; i++) {
        curve2_nodes.push_back(c2.m_residue(i).m_P().coordinates());
    }
    for (i = 0; i < len1; i++) {
        base_positions.push_back(c1.m_residue(i).m_B().coordinates());
    }

    pinang::Vec3d v1(0,0,0);
    pinang::Vec3d v2(0,0,0);
    pinang::Vec3d v3(0,0,0);
    pinang::Vec3d s_tmp(0,0,0);
    pinang::Vec3d t_tmp(0,0,0);

    pinang::Vec3d p_i(0,0,0);   // used in the algorithm of spline fitting;
    pinang::Vec3d b_r(0,0,0);   // used in the algorithm of spline fitting;
    double t1, t2;

    // ==================== curve 1 calculation ====================
    len1 = curve1_nodes.size();
    len2 = curve2_nodes.size();
    for (i = 0; i < len1; i++) {
        if (i == 0)
        {
            s_tmp = curve1_nodes[i+1] - curve1_nodes[i];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve1_tangents.push_back(t_tmp);
        } else if (i == len1-1) {
            s_tmp = curve1_nodes[i] - curve1_nodes[i-1];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve1_tangents.push_back(t_tmp);
        } else {
            v1 = curve1_nodes[i+1] - curve1_nodes[i];
            v2 = curve1_nodes[i] - curve1_nodes[i-1];
            t1 = v1.sqr_norm();
            t2 = v2.sqr_norm();
            s_tmp = v1*t2 + v2*t1;
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve1_tangents.push_back(t_tmp);
        }

        t_tmp = curve2_nodes[len2-1-i] - curve1_nodes[i]; // for cross_vectors;
        cross_vecs.push_back(t_tmp);
    }
    // ~~~~~~~~~~ correction of the first and the last vectors ~~~~~~~~~~
    t1 = (curve1_tangents[1] * curve1_tangents[2]) * 2;
    v1 = curve1_tangents[1] * t1;
    curve1_tangents[0] = v1 - curve1_tangents[2]; // fixed tangent 0;

    t1 = (curve1_tangents[len1-3] * curve1_tangents[len1-2]) * 2;
    v1 = curve1_tangents[len1-2] * t1;
    curve1_tangents[len1-1] = v1 - curve1_tangents[len1-3]; // fixed last tangent;

    // ~~~~~~~~~~ spline fitting ~~~~~~~~~~
    for (i = 0; i < len1 - 1; i++) {
        p_i = curve1_nodes[i];
        t_tmp = curve1_nodes[i+1] - curve1_nodes[i];
        d = t_tmp.norm();
        v1 = curve1_tangents[i];
        v2 = curve1_tangents[i+1];
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r ) + v2 * d * (r*r*r - r*r);
            curve1_dots.push_back(b_r);
        }
    }
    p_i = curve1_nodes[len1-1];
    curve1_dots.push_back(p_i);

    // ==================== curve 2 calculation ====================
    for (i = 0; i < len2; i++) {
        if (i == 0)
        {
            s_tmp = curve2_nodes[i+1] - curve2_nodes[i];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve2_tangents.push_back(t_tmp);
        } else if (i == len2-1) {
            s_tmp = curve2_nodes[i] - curve2_nodes[i-1];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve2_tangents.push_back(t_tmp);
        } else {
            v1 = curve2_nodes[i+1] - curve2_nodes[i];
            v2 = curve2_nodes[i] - curve2_nodes[i-1];
            t1 = v1.sqr_norm();
            t2 = v2.sqr_norm();
            s_tmp = v1*t2 + v2*t1;
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve2_tangents.push_back(t_tmp);
        }
    }

    // ~~~~~~~~~~ correction of the first and the last vectors ~~~~~~~~~~
    t1 = (curve2_tangents[1] * curve2_tangents[2]) * 2;
    v1 = curve2_tangents[1] * t1;
    curve2_tangents[0] = v1 - curve2_tangents[2]; // fixed tangent 0;

    t1 = (curve2_tangents[len1-3] * curve2_tangents[len1-2]) * 2;
    v1 = curve2_tangents[len1-2] * t1;
    curve2_tangents[len1-1] = v1 - curve2_tangents[len1-3]; // fixed last tangent;

    // ~~~~~~~~~~ spline fitting ~~~~~~~~~~
    for (i = 0; i < len2 - 1; i++) {
        p_i = curve2_nodes[i];
        t_tmp = curve2_nodes[i+1] - curve2_nodes[i];
        d = t_tmp.norm();
        v1 = curve2_tangents[i];
        v2 = curve2_tangents[i+1];
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r ) + v2 * d * (r*r*r - r*r);
            curve2_dots.push_back(b_r);
        }
    }
    p_i = curve2_nodes[len2-1];
    curve2_dots.push_back(p_i);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ==================== axis calculation ====================

    pinang::Vec3d H(0,0,0);     // for use in the algorithm of axis calculation;
    pinang::Vec3d P1(0,0,0);    // for use in the algorithm of axis calculation;
    pinang::Vec3d P2(0,0,0);    // for use in the algorithm of axis calculation;
    pinang::Vec3d n1(0,0,0);    // for use in the algorithm of axis calculation;
    pinang::Vec3d n2(0,0,0);    // for use in the algorithm of axis calculation;
    for (i = 0; i < len1; i++) {
        if (i == len1-1) {
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
        r = 2.12 * v1.norm() / p_D_phosphate;
        t_tmp = P1 + (n1 * r);
        axis_nodes.push_back(t_tmp);

        // Calculation of base rise!
        P1 = base_positions[i];
        P2 = base_positions[i+1];
        d = (P2 - P1) * H;
        base_rise.push_back(d);

        // r = helix width!  Also empirical estimation!
        r = 23.0 * v1.norm() / p_D_phosphate;
        helix_width.push_back(r);

        // -------------------- Groove Width! --------------------
        int m = 0, n = 0;
        double major_w = 1000;
        double minor_w = 1000;
        // ---------- major groove ----------
        for (m = i*10 - 50; m < i*10; m++) {
            int proj = 1000;
            double t_dist = 1000;
            if (m < 0)
                continue;
            else {
                v1 = curve1_dots[m] - t_tmp; // vector from P[m] atom to axis center[i];
                n2 = v1 ^ n1;                 // normal vector of the plane;
                n2 = n2 * (1/n2.norm());
                for (n = (len2-i-5)*10; n < (len2-i)*10; n++) {
                    if (n < 0)
                        continue;
                    else {
                        // vector from P atom from the other strand to axis center;
                        v2 = curve2_dots[n] - t_tmp;
                        v2 = v2*(1/v2.norm());
                        if (abs(n2 * v2) <= 0.05) // almost in the same plane!
                        {
                            proj = 0;
                            break;
                        }
                    }
                }
                if (proj == 0)
                    t_dist = vec_distance(curve1_dots[m], curve2_dots[n]);
                if (t_dist < major_w)
                    major_w = t_dist;
            }
        }
        // std::cout << i << "  " << major_w << std::endl;
        major_groove_width.push_back(major_w);

        // ---------- minor groove ----------
        for (m = i * 10 ; m < i * 10 + 60; m++) {
            int proj = 1000;
            double t_dist = 1000;
            if (m > len1 * 10 - 1)
                continue;
            else {
                v1 = curve1_dots[m] - t_tmp; // vector from P[m] atom to axis center[i];
                n2 = v1 ^ n1;                 // normal vector of the plane;
                n2 = n2 * (1/n2.norm());
                for (n = (len2-i) * 10; n < (len2 - i + 6) * 10; n++) {
                    if (n  > len2*10 -1)
                        continue;
                    else {
                        // vector from P atom from the other strand to axis center;
                        v2 = curve2_dots[n] - t_tmp;
                        v2 = v2*(1/v2.norm());
                        if (abs(n2 * v2) <= 0.05) // almost in the same plane!
                        {
                            proj = 0;
                            break;
                        }
                    }
                }
                if (proj == 0)
                    t_dist = vec_distance(curve1_dots[m], curve2_dots[n]);
                if (t_dist < minor_w)
                    minor_w = t_dist;
            }
        }
        // std::cout << i << "  " << minor_w << std::endl;
        minor_groove_width.push_back(minor_w);

    }

    // axis spline fitting!
    for (i = 0; i < axis_nodes.size() - 1; i++) {
        p_i = axis_nodes[i];
        t_tmp = axis_nodes[i+1] - axis_nodes[i];
        d = t_tmp.norm();
        v1 = axis_directions[i];
        v2 = axis_directions[i+1];
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r ) + v2 * d * (r*r*r - r*r);
            axis_dots.push_back(b_r);
        }
    }
    p_i = axis_nodes[axis_nodes.size()-1];
    axis_dots.push_back(p_i);

    // ~~~~~~~~~~~~~~~~~~~~ calculating base per turn ~~~~~~~~~~~~~~~~~~~~
    for (i = 0; i < curve1_tangents.size()-10; i++) {
        t_tmp = curve1_tangents[i];
        double angle = pinang::g_pi;
        int best_fit = 0;
        for (int j = i * 10 + 90; j < i * 10 + 110; j++) {
            if (j > curve1_dots.size()-2)
                break;
            v1 = curve1_dots[j+1] - curve1_dots[j-1];
            s_tmp = v1 * (1/v1.norm());
            double ang = abs(vec_angle(t_tmp, s_tmp));
            if (ang < angle)
            {
                angle = ang;
                best_fit = j;
            }
        }
        bases_per_turn.push_back(best_fit/10.0-i);
    }


    // ============================ Output to PDB ============================
    int k = curve1_dots.size();
    int l = curve1_nodes.size();
    for (i = 0; i < curve1_dots.size(); i++) {
        back_file << std::setw(6) << "HETATM"
                   << std::setw(5) << i+1 << " "
                   << std::setw(4) << "O   "
                   << std::setw(1) << " "
                   << std::setw(3) << "CUR" << " "
                   << std::setw(1) << "A"
                   << std::setw(4) << i/10+2
                   << std::setw(1) << " " << "   "
                   << curve1_dots[i]
                   << std::endl;
    }
    for (i = 0; i < curve2_dots.size(); i++) {
        back_file << std::setw(6) << "HETATM"
                   << std::setw(5) << i+1+k << " "
                   << std::setw(4) << "O   "
                   << std::setw(1) << " "
                   << std::setw(3) << "CUR" << " "
                   << std::setw(1) << "A"
                   << std::setw(4) << i/10+l+3
                   << std::setw(1) << " " << "   "
                   << curve2_dots[i]
                   << std::endl;
    }
    // connecting points!
    for (i = 0; i < curve1_dots.size()-1; i++) {
        back_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1
                   << std::setw(5) << i+2
                   << std::endl;
    }
    for (i = 0; i < curve2_dots.size()-1; i++) {
        back_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1+k
                   << std::setw(5) << i+2+k
                   << std::endl;
    }


    // ~~~~~~~~~~~~~~~~~~~~ output axis pdb ~~~~~~~~~~~~~~~~~~~~
    for (i = 0; i < axis_dots.size(); i++) {
        axis_file << std::setw(6) << "HETATM"
                   << std::setw(5) << i+1 << " "
                   << std::setw(4) << "A   "
                   << std::setw(1) << " "
                   << std::setw(3) << "AXS" << " "
                   << std::setw(1) << "A"
                   << std::setw(4) << i/10+1
                   << std::setw(1) << " " << "   "
                   << axis_dots[i]
                   << std::endl;
    }
    for (i = 0; i < axis_dots.size()-1; i++) {
        axis_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1
                   << std::setw(5) << i+2
                   << std::endl;
    }

    // =================== output base rise, helix width... ===================
    out_file << "# Base rise and helix width:" << std::endl;
    out_file << "#    i - j  " << std::setw(8) << "b_rise"
             << "  " << setw(8) << "h_width" << std::endl;
    for (i = 0; i < len1; i++) {
        out_file << std::setw(6) << i+1
                 << std::setw(4) << i+2 << "  "
                 << std::setw(8)
                 << std::setiosflags(std::ios_base::fixed)
                 << std::setprecision(2)
                 << base_rise[i] << "  "
                 << std::setw(8)
                 << helix_width[i]
                 << std::endl;
    }
    out_file << std::endl << "# Bases-per-turn" << std::endl;
    out_file << std::setw(6) << "# turn " << std::setw(8) << "bpt"
             << std::endl;
    for (i = 0; i < bases_per_turn.size(); i++) {
        out_file << std::setw(6) << i+1 << " "
                 << std::setw(8)
                 << std::setiosflags(std::ios_base::fixed)
                 << std::setprecision(2)
                 << bases_per_turn[i]
                 << std::endl;
    }
    out_file << std::endl << "# Major groove width:" << std::endl
             << std::setw(6) << "#    i" << std::setw(8) << "width"
             << std::endl;
    for (i = 3; i < major_groove_width.size()-3; i++) {
        if (major_groove_width[i] > 25)
            continue;
        out_file << std::setw(6) << i+1
                 << std::setw(8)
                 << std::setiosflags(std::ios_base::fixed)
                 << std::setprecision(2)
                 << major_groove_width[i]
                 << std::endl;
    }
    out_file << std::endl << "# Minor groove width:" << std::endl
             << std::setw(6) << "#    i" << std::setw(8) << "width"
             << std::endl;
    for (i = 3; i < minor_groove_width.size()-3; i++) {
        if (minor_groove_width[i] > 18)
            continue;
        out_file << std::setw(6) << i+1
                 << std::setw(8)
                 << std::setiosflags(std::ios_base::fixed)
                 << std::setprecision(2)
                 << minor_groove_width[i]
                 << std::endl;
    }


    back_file.close();
    axis_file.close();
    out_file.close();

    return 0;
}

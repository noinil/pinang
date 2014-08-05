#include "PDB.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
    std::string infilename = "../data/dna_cg.pdb";
    pinang::PDB pdb1(infilename);

    std::string curve_name = "_curve.pdb";
    std::ofstream curve_file(curve_name.c_str());

    std::string axis_name = "_axis.pdb";
    std::ofstream axis_file(axis_name.c_str());
    std::ofstream test_file("test.pdb");

    std::vector<pinang::Vec3d> curve1_dots;  // Interpolation points;
    std::vector<pinang::Vec3d> curve1_nodes; // "P" atoms in the 1st backbone;
    std::vector<pinang::Vec3d> curve1_tangents; // Tangents at "P" atoms;

    std::vector<pinang::Vec3d> curve2_dots;
    std::vector<pinang::Vec3d> curve2_nodes; // "P " atoms in the 2nd backbone;
    std::vector<pinang::Vec3d> curve2_tangents;

    std::vector<pinang::Vec3d> cross_vecs;

    std::vector<pinang::Vec3d> axis_dots;
    std::vector<pinang::Vec3d> axis_nodes;
    std::vector<pinang::Vec3d> axis_directions; // Axis direction vectors;

    std::vector<double> base_rise;
    std::vector<double> helix_width;
    double d = 0;               // tmp for base_rise;
    double r = 0;               // tmp for helix width;

    int i  = 0;
    int len1 = pdb1.m_model(0).m_chain(0).m_chain_length();
    int len2 = pdb1.m_model(0).m_chain(1).m_chain_length();
    for (i = 1; i < len1; i++) {
        curve1_nodes.push_back(pdb1.m_model(0).m_chain(0).m_residue(i).m_P().coordinates());
    }
    for (i = 1; i < len2; i++) {
        curve2_nodes.push_back(pdb1.m_model(0).m_chain(1).m_residue(i).m_P().coordinates());
    }


    pinang::Vec3d v1(0,0,0);
    pinang::Vec3d v2(0,0,0);
    pinang::Vec3d v3(0,0,0);

    // ==================== curve 1 calculation ====================
    pinang::Vec3d s_tmp(0,0,0);
    pinang::Vec3d t_tmp(0,0,0);
    double t1, t2;
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

        t_tmp = curve2_nodes[len2-1-i] - curve1_nodes[i];
        cross_vecs.push_back(t_tmp);
    }

    pinang::Vec3d p_i(0,0,0);
    pinang::Vec3d b_r(0,0,0);
    for (i = 0; i < len1 - 1; i++) {
        // std::cout << i << " -- nodes:  "
        //           << curve1_nodes[i] << " -- tangents: "
        //           << curve1_tangents[i] << std::endl;
        p_i = curve1_nodes[i];
        t_tmp = curve1_nodes[i+1] - curve1_nodes[i];
        d = t_tmp.norm();
        v1 = curve1_tangents[i];
        v2 = curve1_tangents[i+1];
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r )
                + v2 * d * (r*r*r - r*r);
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

    for (i = 0; i < len2 - 1; i++) {
        // std::cout << i << " -- nodes:  "
        //           << curve2_nodes[i] << " -- tangents: "
        //           << curve2_tangents[i] << std::endl;
        p_i = curve2_nodes[i];
        t_tmp = curve2_nodes[i+1] - curve2_nodes[i];
        d = t_tmp.norm();
        v1 = curve2_tangents[i];
        v2 = curve2_tangents[i+1];
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r )
                + v2 * d * (r*r*r - r*r);
            curve2_dots.push_back(b_r);
        }
    }
    p_i = curve2_nodes[len2-1];
    curve2_dots.push_back(p_i);

    // ==================== axis calculation ====================

    pinang::Vec3d H(0,0,0);
    pinang::Vec3d P1(0,0,0);
    pinang::Vec3d P2(0,0,0);
    pinang::Vec3d n1(0,0,0);
    pinang::Vec3d n2(0,0,0);
    for (i = 0; i < len1; i++) {
        if (i == 0)
        {
            v1 = cross_vecs[i];
            v2 = cross_vecs[i+1];
            t_tmp = v1^v2;
            H = t_tmp * (1/t_tmp.norm());
            P1 = curve1_nodes[i] + (v1 * 0.5);
            P2 = curve1_nodes[i+1] + (v2 * 0.5);
            t_tmp = H ^ v1;
            n1 = t_tmp * (1/t_tmp.norm());
            t_tmp = H ^ v2;
            n2 = t_tmp * (1/t_tmp.norm());
        } else if (i == len1-1) {
            v2 = cross_vecs[i-1];
            v1 = cross_vecs[i];
            t_tmp = v1^v2;
            H = t_tmp * (1/t_tmp.norm());
            P2 = curve1_nodes[i-1] + (v2 * 0.5);
            P1 = curve1_nodes[i] + (v1 * 0.5);
            t_tmp = v1 ^ H;
            n1 = t_tmp * (1/t_tmp.norm());
            t_tmp = v2 ^ H;
            n2 = t_tmp * (1/t_tmp.norm());
        } else {
            v1 = cross_vecs[i];
            v2 = cross_vecs[i+1];
            t_tmp = v1^v2;
            H = t_tmp * (1/t_tmp.norm());
            P1 = curve1_nodes[i] + (v1 * 0.5);
            P2 = curve1_nodes[i+1] + (v2 * 0.5);
            t_tmp = H ^ v1;
            n1 = t_tmp * (1/t_tmp.norm());
            t_tmp = H ^ v2;
            n2 = t_tmp * (1/t_tmp.norm());
        }
        axis_directions.push_back(H);
        s_tmp = P2 - P1;
        d = s_tmp * H;
        // std::cout << " Base rise at " << i << " : " << d << std::endl;
        base_rise.push_back(d);
        // r = ( d*d - s_tmp.sqr_norm())/(2 * (s_tmp * n2));
        r = 2.3 * v1.norm() / 17.45;
        std::cout << " Node " << i+1
                  << " hw " << cross_vecs[i].norm()
                  << " P-dis: " << s_tmp.norm()
                  << " base rise: " << d
                  << std::endl;
        helix_width.push_back(r);
        t_tmp = P1 + (n1 * r);
        axis_nodes.push_back(t_tmp);
    }

    for (i = 0; i < axis_nodes.size() - 1; i++) {
        // std::cout << i << " -- nodes:  "
        //           << axis_nodes[i] << " -- tangents: "
        //           << axis_tangents[i] << std::endl;
        p_i = axis_nodes[i];
        t_tmp = axis_nodes[i+1] - axis_nodes[i];
        d = t_tmp.norm();
        v1 = axis_directions[i];
        v2 = axis_directions[i+1];
        for (int j = 0; j < 10; j++) {
            r = j * 0.1;
            b_r = p_i + v1 * d * ( r - 2*r*r + r*r*r)
                + t_tmp * ( 3*r*r - 2*r*r*r )
                + v2 * d * (r*r*r - r*r);
            axis_dots.push_back(b_r);
        }
    }
    p_i = axis_nodes[axis_nodes.size()-1];
    axis_dots.push_back(p_i);

    // ============================ Output to PDB ============================
    for (i = 0; i < curve1_dots.size(); i++) {
        curve_file << std::setw(6) << "HETATM"
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
    int k = curve1_dots.size();
    int l = curve1_nodes.size();
    for (i = 0; i < curve2_dots.size(); i++) {
        curve_file << std::setw(6) << "HETATM"
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

    for (i = 0; i < curve1_dots.size()-1; i++) {
        curve_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1
                   << std::setw(5) << i+2
                   << std::endl;
    }
    for (i = 0; i < curve2_dots.size()-1; i++) {
        curve_file << std::setw(6) << "CONECT"
                   << std::setw(5) << i+1+k
                   << std::setw(5) << i+2+k
                   << std::endl;
    }


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

    for (i = 0; i < axis_nodes.size(); i++) {
        test_file << std::setw(6) << "HETATM"
                  << std::setw(5) << i+1 << " "
                  << std::setw(4) << "A   "
                  << std::setw(1) << " "
                  << std::setw(3) << "AXS" << " "
                  << std::setw(1) << "A"
                  << std::setw(4) << i/10+1
                  << std::setw(1) << " " << "   "
                  << axis_nodes[i]
                  << std::endl;
    }
    for (i = 0; i < axis_nodes.size()-1; i++) {
        test_file << std::setw(6) << "CONECT"
                  << std::setw(5) << i+1
                  << std::setw(5) << i+2
                  << std::endl;
    }

    curve_file.close();
    axis_file.close();

    return 0;
}

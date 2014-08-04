#include "PDB.h"
#include "vec3d.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
    std::string infilename = "../data/dna0.pdb";
    pinang::PDB pdb1(infilename);

    std::string curve_name = "_curve.pdb";
    std::ofstream curve_file(curve_name.c_str());

    std::vector<pinang::Vec3d> curve1_dots;
    std::vector<pinang::Vec3d> curve1_nodes;
    std::vector<pinang::Vec3d> curve1_tangents;

    std::vector<pinang::Vec3d> curve2_dots;
    std::vector<pinang::Vec3d> curve2_nodes;
    std::vector<pinang::Vec3d> curve2_tangents;

    int i  = 0;
    for (i = 1; i < pdb1.m_model(0).m_chain(0).m_chain_length(); i++) {
        curve1_nodes.push_back(pdb1.m_model(0).m_chain(0).m_residue(i).m_P().coordinates());
    }
    for (i = 1; i < pdb1.m_model(0).m_chain(1).m_chain_length(); i++) {
        curve2_nodes.push_back(pdb1.m_model(0).m_chain(1).m_residue(i).m_P().coordinates());
    }

    // ==================== curve 1 calculation ====================
    pinang::Vec3d s_tmp(0,0,0);
    pinang::Vec3d t_tmp(0,0,0);
    pinang::Vec3d v1(0,0,0);
    pinang::Vec3d v2(0,0,0);
    double t1, t2;
    for (i = 0; i < curve1_nodes.size(); i++) {
        if (i == 0)
        {
            s_tmp = curve1_nodes[i+1] - curve1_nodes[i];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve1_tangents.push_back(t_tmp);
        } else if (i == curve1_nodes.size()-1) {
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
    }

    double r = 0;
    double d = 0;
    pinang::Vec3d p_i(0,0,0);
    pinang::Vec3d b_r(0,0,0);
    for (i = 0; i < curve1_nodes.size() - 1; i++) {
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
    p_i = curve1_nodes[curve1_nodes.size()-1];
    curve1_dots.push_back(p_i);

    // ==================== curve 2 calculation ====================
    for (i = 0; i < curve2_nodes.size(); i++) {
        if (i == 0)
        {
            s_tmp = curve2_nodes[i+1] - curve2_nodes[i];
            t_tmp = s_tmp * (1 / s_tmp.norm());
            curve2_tangents.push_back(t_tmp);
        } else if (i == curve2_nodes.size()-1) {
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

    for (i = 0; i < curve2_nodes.size() - 1; i++) {
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
    p_i = curve2_nodes[curve2_nodes.size()-1];
    curve2_dots.push_back(p_i);

    // ============================ Output to PDB ============================
    for (i = 0; i < curve1_dots.size(); i++) {
        curve_file << std::setw(6) << "HETATM"
                   << std::setw(5) << i+1 << " "
                   << std::setw(4) << "O   "
                   << std::setw(1) << " "
                   << std::setw(3) << "CUR" << " "
                   << std::setw(1) << "A"
                   << std::setw(4) << i/10+1
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
                   << std::setw(4) << i/10+l+2
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

}

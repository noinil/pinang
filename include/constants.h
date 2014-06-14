// -*-c++-*-

#ifndef PINANG_CONSTANTS_H_
#define PINANG_CONSTANTS_H_

namespace pinang {

    double g_cutoff = 6.5;

    const long double g_pi = 3.14159265358979323846;

    // Units.
    const double u_mass = 1.0;

    // CG MD Parameters.
    const double p_K_bond = 100.0;
    const double p_K_angle = 20.0;
    const double p_K_dihedral_1 = 1.0;
    const double p_K_dihedral_3 = 0.5;
    const double p_K_native = 1;
    const double p_K_nonnative = 1;
}

#endif

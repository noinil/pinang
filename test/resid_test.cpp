#include <iostream>
#include "residue.h"

using namespace std;

int main(int argc, char *argv[])
{
    pinang::Atom a, b, c;
    pinang::Vec3d ca, cb, cc, dd;

    ca = pinang::Vec3d(1, 2, 3);
    cb = pinang::Vec3d(1, 0, 0);
    cc = pinang::Vec3d(0, 1, 0);
    dd = (ca^cb) + cc;

    std::string sid("LOS");
    // double o;
    // int i = 0;

    pinang::Residue r_0;

    a.set_serial(1);
    a.set_atom_name("O1");
    a.set_resid_name("ASP");
    a.set_resid_index(100);
    a.set_coords(ca);
    a.set_segment_ID(sid);

    b.set_serial(2);
    b.set_atom_name("N");
    b.set_resid_name("ASP");
    b.set_resid_index(100);
    b.set_coords(cb);

    c.set_serial(3);
    c.set_atom_name("CA");
    c.set_resid_name("ASP");
    c.set_resid_index(100);
    c.set_coords(cc);

    r_0.set_resid_name("ASP");
    r_0.set_chain_ID('A');
    r_0.set_resid_index(1);

    // std::cout << "the first atom of residue:"
    //           << r_0.m_atom(1).atom_name()
    //           << endl;
    r_0.add_atom(a);
    r_0.add_atom(b);
    r_0.add_atom(c);

    std::cout << "the first atom of residue: "
              << r_0.m_atom(0).atom_name()
              << endl;

    std::cout << "the first atom of residue: "
              << r_0.m_atom(1).atom_name()
              << endl;

    std::cout << "the first atom of residue: "
              << r_0.m_atom(2).atom_name()
              << endl;

    return 0;
}

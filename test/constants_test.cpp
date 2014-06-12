#include <iostream>
#include "constants.h"
#include "vec3d.h"
#include "contact.h"

using namespace std;

int main(int argc, char *argv[])
{
    // pinang::g_cutoff = 1.0;

    // std::cout << "PINANG::CUTOFF = " << pinang::g_cutoff << std::endl;

    pinang::Vec3d a0(0,0,0);
    pinang::Vec3d a1(1,0,0);
    pinang::Vec3d a2(0,1,0);
    pinang::Vec3d a3(0,0,12);

    pinang::Contact c(1,2,vec_distance(a1,a2));

    int t = c.is_contact();

    std::cout << "a1 and a2 are in contact? " << t << std::endl;

    c.set(1, 3, vec_distance(a1, a3));
    t = c.is_contact();

    std::cout << "a1 and a3 are in contact? " << t << std::endl;

    return 0;
}

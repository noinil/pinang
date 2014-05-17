#include <iostream>
#include <cmath>
#include "binang/vec3d.h"

using namespace std;

int main(int argc, char *argv[])
{
    binang::Vec3d a(1,0,0);
    binang::Vec3d b(0,1,0);
    binang::Vec3d c(0,1,1);
    binang::Vec3d d;
    double m;

    std::cout << "a = " << a << endl;
    std::cout << "b = " << b << endl;
    std::cout << "c = " << c << endl;

    d = a+b;
    std::cout << "d = a + b = " << d << endl;

    d = a - b;
    cout << "d = a - b = " << d << endl;

    d = c * 3 + b;
    cout << "d = c * 3 + b = " << d << endl;
    cout << "d.norm() = " << d.norm() << endl;

    double dist = Vec_distance(a,b);
    cout << "distance from a to b is " << dist << endl;
    cout << "distance squared from a to b is " << pow(dist, 2) << endl;

    m = c*b;
    cout << "a * b = " << m << endl;
    d = c^b;
    cout << "c x b = " << d << endl;
    d = b^c;
    cout << "b x c = " << d << endl;

    return 0;
}

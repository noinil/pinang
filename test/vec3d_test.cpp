#include <iostream>
#include <cmath>
#include "vec3d.h"
#include <limits>


using namespace std;

int main(int argc, char *argv[])
{
    pinang::Vec3d a(1,0,0);
    pinang::Vec3d b(0,1,0);
    pinang::Vec3d c(0,1,1);
    pinang::Vec3d d;
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

    double dist = vec_distance(a,b);
    cout << "distance from a to b is " << dist << endl;
    cout << "distance squared from a to b is " << pow(dist, 2) << endl;

    m = c*b;
    cout << "a * b = " << m << endl;
    d = c^b;
    cout << "c x b = " << d << endl;
    d = b^c;
    cout << "b x c = " << d << endl;

    cout << endl << "test operator[]: a[0] = " << a[0] << endl;
    cout << endl << "angle between a and b = " << vec_angle_deg(a, b);
    cout << endl << "angle between c and b = " << vec_angle_deg(b, c);
    cout << endl << "angle between a and c = " << vec_angle_deg(a, c);

    std::cout << endl
              << std::setprecision(std::numeric_limits<long double>::digits10)
              << std::setw(20) << pinang::g_pi << std::endl;
    return 0;
}

#include <iostream>
#include "binang/atom.h"
 #include "binang/vec3d.h"

using namespace std;

int main(int argc, char *argv[])
{
    binang::Atom a, b, c;
    binang::Vec3d ca, cb, cc, dd;

    ca = binang::Vec3d(1, 2, 3);
    cb = binang::Vec3d(1, 0, 0);
    cc = binang::Vec3d(0, 1, 0);

    std::string sid("LOS");
    double o;
    int i = 0;

    cout << "serials:" << endl;
    cout << "a: " << a.serial() << endl;
    cout << "b: " << b.serial() << endl;
    cout << "c: " << c.serial() << endl;

    a.set_serial(i+1);
    b.set_serial(2);
    c.set_serial(3);

    cout << "serials:" << endl;
    cout << "a: " << a.serial() << endl;
    cout << "b: " << b.serial() << endl;
    cout << "c: " << c.serial() << endl;

    o = 0.6;
    a.set_occupancy(o);
    b.set_occupancy(.8);
    c.set_occupancy(1);

    cout << "occupancy:" << endl;
    cout << "a: " << a.occupancy() << endl;
    cout << "b: " << b.occupancy() << endl;
    cout << "c: " << c.occupancy() << endl;

    cout << "seg_ID:" << endl;
    cout << "a: " << a.segment_ID() << endl;
    cout << "b: " << b.segment_ID() << endl;
    cout << "c: " << c.segment_ID() << endl;
    c.segment_ID() = "bbb";
    cout << "c: " << c.segment_ID() << endl;

    a.set_segment_ID(sid);
    b.set_segment_ID("TST");
    c.set_segment_ID("aaa");

    cout << "seg_ID:" << endl;
    cout << "a: " << a.segment_ID() << endl;
    cout << "b: " << b.segment_ID() << endl;
    cout << "c: " << c.segment_ID() << endl;
    c.segment_ID() = "bbb";
    cout << "c: " << c.segment_ID() << endl;

    cout << "coordinates:" << endl;
    cout << "a: " << a.coordinates() << endl;
    cout << "b: " << b.coordinates() << endl;
    cout << "c: " << c.coordinates() << endl;

    a.set_coords(ca);
    b.set_coords(cb);
    c.set_coords(cc);
    cout << "coordinates:" << endl;
    cout << "a: " << a.coordinates() << endl;
    cout << "b: " << b.coordinates() << endl;
    cout << "c: " << c.coordinates() << endl;

    dd = c.coordinates();
    cout << "d: " << dd << endl;
    dd = ca;
    cout << "d: " << dd << endl;
    cout << "c: " << c.coordinates() << endl;

    // a = b;
    // cout << "a: " << a.coordinates() << endl;
    // cout << "b: " << b.coordinates() << endl;

    return 0;
}

// -*-c++-*-

#ifndef PINANG_CONFORMATION_H
#define PINANG_CONFORMATION_H

#include "constants.h"
#include "vec3d.h"

#include <vector>


namespace pinang {
    class Conformation
    {
    public:
        Conformation();
        Conformation(std::vector<Vec3d> v);

        inline void reset();

        int m_size() {return _n_atom;}

        int set_conformation(std::vector<Vec3d> v);
        inline const Vec3d& atom(int n) const;

        virtual ~Conformation() {_coordinates.clear();}
    protected:
        std::vector<Vec3d> _coordinates;
        int _n_atom;
    };

    Conformation::Conformation()
    {
        _n_atom = 0;
        _coordinates.clear();
    }

    Conformation::Conformation(std::vector<Vec3d> v)
    {
        _coordinates = v;
        _n_atom = _coordinates.size();
    }

    inline void Conformation::reset()
    {
        _n_atom = 0;
        _coordinates.clear();
    }

    int Conformation::set_conformation(std::vector<Vec3d> v)
    {
        int m = v.size();
        if (m != _n_atom && _n_atom > 0)
        {
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~           PINANG :: CONFORMATION           ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << " ERROR: Wrong atom number when set conformation. " << std::endl;
            return 1;
        } else {
            _coordinates = v;
            _n_atom = m;
            return 0;
        }
        // _coordinates = v;
        // _n_atom = m;
        // return 0;
    }

    inline const Vec3d& Conformation::atom(int n) const
    {
        if ( n >= _n_atom || n < 0)
        {
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~           PINANG :: CONFORMATION           ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << " ERROR: Atom index out of range in Conformation. " << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            return _coordinates[n];
        }
    }

}
#endif

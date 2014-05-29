// -*-c++-*-

#ifndef PINANG_RESIDUE_H_
#define PINANG_RESIDUE_H_

#include <iostream>
#include <vector>
#include <cstdlib>

#include "vec3d.h"
#include "atom.h"

namespace pinang {
    class Residue
    {
    public:
        Residue();
        virtual ~Residue() {_atoms.clear();}

        inline void reset();

        inline std::string resid_name() const;
        inline void set_resid_name(const std::string& s);

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline unsigned int resid_index() const;
        inline void set_resid_index(unsigned int i);

        inline Atom& m_atom(unsigned int n);
        inline int add_atom(const Atom& a);

        inline int m_residue_size() const;

        inline Atom& m_C_alpha();

    protected:
        std::string _resid_name;
        char _chain_ID;
        unsigned int _resid_index;
        std::vector<Atom> _atoms;
        int _n_atom;

        Atom _C_alpha;
    };

    /*                _     _
    //  _ __ ___  ___(_) __| |  _ __   __ _ _ __ ___   ___
    // | '__/ _ \/ __| |/ _` | | '_ \ / _` | '_ ` _ \ / _ \
    // | | |  __/\__ \ | (_| | | | | | (_| | | | | | |  __/
    // |_|  \___||___/_|\__,_| |_| |_|\__,_|_| |_| |_|\___|
    */
    inline std::string Residue::resid_name() const
    {
        return _resid_name;
    }
    inline void Residue::set_resid_name(const std::string& s)
    {
        _resid_name = s;
    }


    /*      _           _         ___ ____
    //  ___| |__   __ _(_)_ __   |_ _|  _ \
    // / __| '_ \ / _` | | '_ \   | || | | |
    //| (__| | | | (_| | | | | |  | || |_| |
    // \___|_| |_|\__,_|_|_| |_| |___|____/
    */
    inline char Residue::chain_ID() const
    {
        return _chain_ID;
    }
    inline void Residue::set_chain_ID(char a)
    {
        _chain_ID = a;
    }

    /*                _     _   _           _
    //  _ __ ___  ___(_) __| | (_)_ __   __| | _____  __
    // | '__/ _ \/ __| |/ _` | | | '_ \ / _` |/ _ \ \/ /
    // | | |  __/\__ \ | (_| | | | | | | (_| |  __/>  <
    // |_|  \___||___/_|\__,_| |_|_| |_|\__,_|\___/_/\_\
    */
    inline unsigned int Residue::resid_index() const
    {
        return _resid_index;
    }
    inline void Residue::set_resid_index(unsigned int i)
    {
        _resid_index = i;
    }

    /*                     _
    //  _ __ ___      __ _| |_ ___  _ __ ___
    // | '_ ` _ \    / _` | __/ _ \| '_ ` _ \
    // | | | | | |  | (_| | || (_) | | | | | |
    // |_| |_| |_|___\__,_|\__\___/|_| |_| |_|
    //          |_____|
    */
    inline Atom& Residue::m_atom(unsigned int n)
    {
        if (_atoms.empty())
        {
            std::cerr << "ERROR: No Atoms found in Residue: "
                      << _resid_index << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _atoms.size())
            {
                std::cerr << "ERROR: Atom index out of range in Residue: "
                          << _resid_index << std::endl;
                exit(EXIT_SUCCESS);
            } else {
                return _atoms[n];
            }
        }
    }
    inline int Residue::add_atom(const Atom& a)
    {
        if (a.resid_index() != _resid_index)
        {
            return 1;
        }
        _atoms.push_back(a);
        if (a.atom_name() == "CA")
        {
            _C_alpha = a;
        }
        _n_atom++;
        return 0;
    }

    inline Atom& Residue::m_C_alpha()
    {
        return _C_alpha;
    }

    inline int Residue::m_residue_size() const
    {
        return _n_atom;
    }

    // Residue------------------------------------------------------------------
    inline Residue::Residue()
    {
        _resid_name = "";
        _chain_ID = -1;
        _resid_index = -1;
        _atoms.clear();
        _n_atom = 0;
        _C_alpha.reset();
    }

    inline void Residue::reset()
    {
        _resid_name = "";
        _chain_ID = -1;
        _resid_index = -1;
        _atoms.clear();
        _n_atom = 0;
        _C_alpha.reset();
    }

    inline std::ostream& operator<<(std::ostream& o, Residue& r)
    {
        o << "Residue "
          << std::setw(4) << r.resid_index() << ":  "
          << std::setw(3) << r.resid_name() << std::endl;
        int i = 0;
        for (i = 0; i < r.m_residue_size(); i++) {
            o << r.m_atom(i) << std::endl;
        }
        return o;
    }

    inline double resid_min_distance (Residue& r1, Residue& r2)
    {
        int i, j;
        double d = -1;           // distance;
        double f = 0;           // tmp distance;
        for (i = 0; i < r1.m_residue_size(); i++) {
            for (j = 0; j < r2.m_residue_size(); j++) {
                f = atom_distance(r1.m_atom(i), r2.m_atom(j));
                if (d < 0 || d > f)
                {
                    d = f;
                }
            }
        }
        return d;
    }

}

#endif

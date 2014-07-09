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
        inline char short_name() const;
        void set_resid_name(const std::string& s);

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline unsigned int resid_index() const;
        inline void set_resid_index(unsigned int i);

        inline double resid_charge() const;
        inline void set_resid_charge(double c);

        inline double resid_mass() const;
        inline void set_resid_mass(double c);

        inline Atom& m_atom(unsigned int n);
        inline int add_atom(const Atom& a);

        inline int m_residue_size() const;

        inline Atom& m_C_alpha();

    protected:
        std::string _resid_name;
        char _short_name;
        char _chain_ID;
        unsigned int _resid_index;
        std::vector<Atom> _atoms;
        int _n_atom;
        double _charge;
        double _mass;

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
    inline char Residue::short_name() const
    {
        return _short_name;
    }
    void Residue::set_resid_name(const std::string& s)
    {
        _resid_name = s;
        char c = s[0];

        switch (c) {
        case 'A':
            if (_resid_name == "ARG") {_short_name = 'R'; _charge = 1.0; _mass = 156.19;}
            else if (_resid_name == "ASP") {_short_name = 'D'; _charge = -1.0; _mass = 115.09;}
            else if (_resid_name == "ASN") {_short_name = 'N'; _mass = 114.11;}
            else if (_resid_name == "ALA") {_short_name = 'A'; _mass = 71.09;}
            else if (_resid_name == "A") _short_name = 'A';
            break;
        case 'C':
            if (_resid_name == "CYS") {_short_name = 'C';  _mass = 103.15;}
            else if (_resid_name == "C") _short_name = 'C';
            else if (_resid_name == "CA") {_short_name = 'c'; _charge = 2.0; _mass = 40.08;}
            break;
        case 'D':
            if (_resid_name == "DA") _short_name = 'A';
            else if (_resid_name == "DC") _short_name = 'C';
            else if (_resid_name == "DG") _short_name = 'G';
            else if (_resid_name == "DT") _short_name = 'T';
            break;
        case 'G':
            if (_resid_name == "GLU") {_short_name = 'E'; _charge = -1.0; _mass = 129.12;}
            else if (_resid_name == "GLN") {_short_name = 'Q';  _mass = 128.14;}
            else if (_resid_name == "GLY") {_short_name = 'G';  _mass = 57.05;}
            else if (_resid_name == "G") _short_name = 'G';
            break;
        case 'H':
            if (_resid_name == "HIS") {_short_name = 'H';  _mass = 137.14; _charge = 1.0;}
            else if (_resid_name == "HOH") _short_name = 'w';
            break;
        case 'I':
            if (_resid_name == "ILE") {_short_name = 'I';  _mass = 113.16;}
            break;
        case 'L':
            if (_resid_name == "LYS") {_short_name = 'K'; _charge = 1.0; _mass = 128.17;}
            else if (_resid_name == "LEU") {_short_name = 'L';  _mass = 113.16;}
            break;
        case 'M':
            if (_resid_name == "MET") {_short_name = 'M';  _mass = 131.19;}
            break;
        case 'P':
            if (_resid_name == "PRO") {_short_name = 'P';  _mass = 97.12;}
            else if (_resid_name == "PHE") {_short_name = 'F';  _mass = 147.18;}
            break;
        case 'R':
            if (_resid_name == "RA") _short_name = 'A';
            else if (_resid_name == "RU") _short_name = 'U';
            else if (_resid_name == "RC") _short_name = 'C';
            else if (_resid_name == "RG") _short_name = 'G';
            break;
        case 'S':
            if (_resid_name == "SER") {_short_name = 'S';  _mass = 87.08;}
            else if (_resid_name == "SEC") {_short_name = 'U';}
            break;
        case 'T':
            if (_resid_name == "TYR") {_short_name = 'Y';  _mass = 163.18;}
            else if (_resid_name == "TRP") {_short_name = 'W'; _mass = 186.21;}
            else if (_resid_name == "THR") {_short_name = 'T'; _mass = 101.11;}
            else if (_resid_name == "T") _short_name = 'T';
            break;
        case 'U':
            if (_resid_name == "U") _short_name = 'U';
            break;
        case 'V':
            if (_resid_name == "VAL") {_short_name = 'V'; _mass = 99.14;}
            break;
        default:
            if (_resid_name == "ZN") {_short_name = 'z'; _charge = 2.0; _mass = 65.37;}
        }

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

    /*                _     _        _
    //  _ __ ___  ___(_) __| |   ___| |__   __ _ _ __ __ _  ___
    // | '__/ _ \/ __| |/ _` |  / __| '_ \ / _` | '__/ _` |/ _ \
    // | | |  __/\__ \ | (_| | | (__| | | | (_| | | | (_| |  __/
    // |_|  \___||___/_|\__,_|  \___|_| |_|\__,_|_|  \__, |\___|
    //                                               |___/
    */
    inline double Residue::resid_charge() const
    {
        return _charge;
    }
    inline void Residue::set_resid_charge(double c)
    {
        _charge = c;
    }

    inline double Residue::resid_mass() const
    {
        return _mass;
    }
    inline void Residue::set_resid_mass(double m)
    {
        _mass = m;
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
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~             PINANG :: Residue              ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << "ERROR: No Atoms found in Residue: "
                      << _resid_index << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n >= _atoms.size())
            {
                std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
                std::cout << " ~             PINANG :: Residue              ~ " << std::endl;
                std::cout << " ============================================== " << std::endl;
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
        if (a.atom_name() == "C3'")
        {
            _C_alpha = a;
        }
        if (a.atom_flag() == "HETATM")
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
        _short_name = '0';
        _chain_ID = -1;
        _resid_index = -1;
        _atoms.clear();
        _n_atom = 0;
        _charge = 0.0;
        _mass = 100.0;

        _C_alpha.reset();
    }

    inline void Residue::reset()
    {
        _resid_name = "";
        _short_name = '0';
        _chain_ID = -1;
        _resid_index = -1;
        _atoms.clear();
        _n_atom = 0;
        _charge = 0.0;
        _mass = 100.0;

        _C_alpha.reset();
    }

    /*              _               __
    //   ___  _   _| |_ ___ _ __   / _|_   _ _ __   ___
    //  / _ \| | | | __/ _ \ '__| | |_| | | | '_ \ / __|
    // | (_) | |_| | ||  __/ |    |  _| |_| | | | | (__
    //  \___/ \__,_|\__\___|_|    |_|  \__,_|_| |_|\___|
    */
    inline std::ostream& operator<<(std::ostream& o, Residue& r)
    {
        // o << "Residue "
        //   << std::setw(4) << r.resid_index() << ":  "
        //   << std::setw(3) << r.resid_name() << std::endl;
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
            if (r1.m_atom(i).element() == "H")
                continue;
            for (j = 0; j < r2.m_residue_size(); j++) {
                if (r2.m_atom(j).element() == "H")
                    continue;
                f = atom_distance(r1.m_atom(i), r2.m_atom(j));
                if (d < 0 || d > f)
                {
                    d = f;
                }
            }
        }
        return d;
    }

    inline double resid_ca_distance (Residue& r1, Residue& r2)
    {
        double d = -1;           // distance;
        d = atom_distance(r1.m_C_alpha(), r2.m_C_alpha());
        return d;
    }

}

#endif

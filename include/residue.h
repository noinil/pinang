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

        inline chain_t chain_type() const;
        inline void set_chain_type(chain_t ct);

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
        inline Atom& m_P();
        inline Atom& m_S();
        inline Atom& m_B();

        void set_cg_na();

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
        Atom _P;
        Atom _S;
        Atom _B;
        chain_t _chain_type;
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
            if (_resid_name == "ARG") {_short_name = 'R'; _charge = 1.0; _mass = 156.19; _chain_type = protein;}
            else if (_resid_name == "ASP") {_short_name = 'D'; _charge = -1.0;
                _mass = 115.09; _chain_type = protein;}
            else if (_resid_name == "ASN") {_short_name = 'N'; _mass = 114.11; _chain_type = protein;}
            else if (_resid_name == "ALA") {_short_name = 'A'; _mass = 71.09; _chain_type = protein;}
            else if (_resid_name == "A") {_short_name = 'A'; _chain_type = na;_mass = 134.12;}
            break;
        case 'C':
            if (_resid_name == "CYS") {_short_name = 'C';  _mass = 103.15;  _chain_type = protein;}
            else if (_resid_name == "C") {_short_name = 'C'; _chain_type = na;  _mass = 110.09;}
            else if (_resid_name == "CA") {_short_name = 'c'; _charge = 2.0; _mass = 40.08; _chain_type = ion;}
            break;
        case 'D':
            if (_resid_name == "DA") {_short_name = 'A'; _chain_type = DNA; _mass = 134.12;}
            else if (_resid_name == "DC") {_short_name = 'C'; _chain_type = DNA;  _mass = 110.09;}
            else if (_resid_name == "DG") {_short_name = 'G'; _chain_type = DNA;  _mass = 150.12;}
            else if (_resid_name == "DT") {_short_name = 'T'; _chain_type = DNA;  _mass = 125.091;}
            break;
        case 'G':
            if (_resid_name == "GLU") {_short_name = 'E'; _charge = -1.0; _mass = 129.12; _chain_type = protein;}
            else if (_resid_name == "GLN") {_short_name = 'Q';  _mass = 128.14;  _chain_type = protein;}
            else if (_resid_name == "GLY") {_short_name = 'G';  _mass = 57.05; _chain_type = protein;}
            else if (_resid_name == "G") {_short_name = 'G'; _chain_type = na;  _mass = 150.12;}
            break;
        case 'H':
            if (_resid_name == "HIS") {_short_name = 'H';  _mass = 137.14; _charge = 1.0; _chain_type = protein;}
            else if (_resid_name == "HOH") {_short_name = 'w';   _chain_type = water;}
            break;
        case 'I':
            if (_resid_name == "ILE") {_short_name = 'I';  _mass = 113.16; _chain_type = protein;}
            break;
        case 'L':
            if (_resid_name == "LYS") {_short_name = 'K'; _charge = 1.0; _mass = 128.17; _chain_type = protein;}
            else if (_resid_name == "LEU") {_short_name = 'L';  _mass = 113.16; _chain_type = protein;}
            break;
        case 'M':
            if (_resid_name == "MET") {_short_name = 'M';  _mass = 131.19; _chain_type = protein;}
            else if (_resid_name == "MG") {_short_name = 'm';  _mass = 24.305; _charge = 2.0; _chain_type = ion;}
            break;
        case 'P':
            if (_resid_name == "PRO") {_short_name = 'P';  _mass = 97.12; _chain_type = protein;}
            else if (_resid_name == "PHE") {_short_name = 'F';  _mass = 147.18; _chain_type = protein;}
            break;
        case 'R':
            if (_resid_name == "RA") {_short_name = 'A'; _chain_type = RNA; _mass = 134.12;}
            else if (_resid_name == "RU") {_short_name = 'U'; _chain_type = RNA; _mass = 111.08;}
            else if (_resid_name == "RC") {_short_name = 'C'; _chain_type = RNA;  _mass = 110.09;}
            else if (_resid_name == "RG") {_short_name = 'G';  _chain_type = RNA; _mass = 150.12;}
            break;
        case 'S':
            if (_resid_name == "SER") {_short_name = 'S';  _mass = 87.08; _chain_type = protein;}
            else if (_resid_name == "SEC") {_short_name = 'U'; _chain_type = protein;}
            break;
        case 'T':
            if (_resid_name == "TYR") {_short_name = 'Y';  _mass = 163.18; _chain_type = protein;}
            else if (_resid_name == "TRP") {_short_name = 'W'; _mass = 186.21; _chain_type = protein;}
            else if (_resid_name == "THR") {_short_name = 'T'; _mass = 101.11; _chain_type = protein;}
            else if (_resid_name == "T") {_short_name = 'T';  _chain_type = DNA; _mass = 125.091;}
            break;
        case 'U':
            if (_resid_name == "U") {_short_name = 'U'; _chain_type = RNA; _mass = 111.08;}
            break;
        case 'V':
            if (_resid_name == "VAL") {_short_name = 'V'; _mass = 99.14; _chain_type = protein;}
            break;
        default:
            if (_resid_name == "ZN") {_short_name = 'z'; _charge = 2.0; _mass = 65.409; _chain_type = ion;}
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

    /*       _           _         _
    //   ___| |__   __ _(_)_ __   | |_ _   _ _ __   ___
    //  / __| '_ \ / _` | | '_ \  | __| | | | '_ \ / _ \
    // | (__| | | | (_| | | | | | | |_| |_| | |_) |  __/
    //  \___|_| |_|\__,_|_|_| |_|  \__|\__, | .__/ \___|
    //                                 |___/|_|
    */
    inline chain_t Residue::chain_type() const
    {
        return _chain_type;
    }
    inline void Residue::set_chain_type(chain_t a)
    {
        _chain_type = a;
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
        if (a.atom_name() == "CA ")
        {
            _C_alpha = a;
        }
        if (a.atom_name() == "C3'" || a.atom_name() == "S  ")
        {
            _S = a;
        }
        if (a.atom_name() == "P  ")
        {
            _P = a;
        }
        if (a.atom_name() == "N1 " || a.atom_name() == "B  ")
        {
            _B = a;
        }
        if (a.atom_flag() == "HETATM" && a.element() != "H")
        {
            _C_alpha = a;
        }
        _n_atom++;
        return 0;
    }

    /*   ____         _       _
    //  / ___|   __ _| |_ __ | |__   __ _
    // | |      / _` | | '_ \| '_ \ / _` |
    // | |___  | (_| | | |_) | | | | (_| |
    //  \____|  \__,_|_| .__/|_| |_|\__,_|
    //                 |_|
    */
    inline Atom& Residue::m_C_alpha()
    {
        return _C_alpha;
    }

    inline Atom& Residue::m_P()
    {
        return _P;
    }

    inline Atom& Residue::m_S()
    {
        return _S;
    }

    inline Atom& Residue::m_B()
    {
        return _B;
    }

    void Residue::set_cg_na()
    {
        int i = 0;
        Vec3d coor_P(0,0,0);
        Vec3d coor_C5p(0,0,0),
            coor_C4p(0,0,0),
            coor_O4p(0,0,0),
            coor_C1p(0,0,0),
            coor_C2p(0,0,0),
            coor_C3p(0,0,0),
            coor_O2p(0,0,0);
        Vec3d coor_N1(0,0,0),
            coor_C2(0,0,0),
            coor_N3(0,0,0),
            coor_C4(0,0,0),
            coor_C5(0,0,0),
            coor_C6(0,0,0),
            coor_O2(0,0,0),
            coor_N4(0,0,0),
            coor_O4(0,0,0),
            coor_N6(0,0,0),
            coor_O6(0,0,0),
            coor_N2(0,0,0),
            coor_N7(0,0,0),
            coor_C8(0,0,0),
            coor_N9(0,0,0);
        Vec3d com_P(0,0,0);
        Vec3d com_S(0,0,0);
        Vec3d com_B(0,0,0);
        double mass_C = 12.011;
        double mass_O = 15.999;
        double mass_N = 14.001;
        int n_cs=0;
        int n_cb=0;
        int n_os=0;
        int n_ob=0;
        int n_nb=0;

        if (_chain_type != DNA && _chain_type != RNA && _chain_type != na)
        {
            return;
        }
        for (i = 0; i < _n_atom; i++) {
            std::string aname = _atoms[i].atom_name();
            char c = aname[0];
            switch (c) {
            case 'C':
                if (aname == "C5'") {coor_C5p = _atoms[i].coordinates(); n_cs++;}
                else if (aname == "C1'") {coor_C1p = _atoms[i].coordinates(); n_cs++;}
                else if (aname == "C2'") {coor_C2p = _atoms[i].coordinates(); n_cs++;}
                else if (aname == "C3'") {coor_C3p = _atoms[i].coordinates(); n_cs++;}
                else if (aname == "C4'") {coor_C4p = _atoms[i].coordinates(); n_cs++;}
                else if (aname == "C2 ") { coor_C2 = _atoms[i].coordinates(); n_cb++;}
                else if (aname == "C4 ") { coor_C4 = _atoms[i].coordinates(); n_cb++;}
                else if (aname == "C5 ") { coor_C5 = _atoms[i].coordinates(); n_cb++;}
                else if (aname == "C6 ") { coor_C6 = _atoms[i].coordinates(); n_cb++;}
                else if (aname == "C8 ") { coor_C8 = _atoms[i].coordinates(); n_cb++;}
                break;
            case 'O':
                if (aname == "O4'") {coor_O4p = _atoms[i].coordinates(); n_os++;}
                else if (aname == "O2'") {coor_O2p = _atoms[i].coordinates(); n_os++;}
                else if (aname == "O2 ") {coor_O2 = _atoms[i].coordinates(); n_ob++;}
                else if (aname == "O4 ") {coor_O4 = _atoms[i].coordinates(); n_ob++;}
                else if (aname == "O6 ") {coor_O6 = _atoms[i].coordinates(); n_ob++;}
                break;
            case 'N':
                if (aname == "N1 ") {coor_N1 = _atoms[i].coordinates(); n_nb++;}
                else if (aname == "N2 ") {coor_N2 = _atoms[i].coordinates(); n_nb++;}
                else if (aname == "N3 ") {coor_N3 = _atoms[i].coordinates(); n_nb++;}
                else if (aname == "N4 ") {coor_N4 = _atoms[i].coordinates(); n_nb++;}
                else if (aname == "N6 ") {coor_N6 = _atoms[i].coordinates(); n_nb++;}
                else if (aname == "N7 ") {coor_N7 = _atoms[i].coordinates(); n_nb++;}
                else if (aname == "N9 ") {coor_N9 = _atoms[i].coordinates(); n_nb++;}
                break;
            default:
                if (aname == "P  ") coor_P = _atoms[i].coordinates();
            }
        }
        com_P = coor_P;
        com_S = ( (coor_C1p + coor_C2p + coor_C3p + coor_C4p + coor_C5p)
                  * mass_C
                  + ( coor_O2p + coor_O4p ) * mass_O ) * (1/( n_os * mass_O + n_cs * mass_C ));
        com_B = ( (coor_N1 + coor_N2 + coor_N3 + coor_N4 + coor_N6 + coor_N7 + coor_N9)
                  * mass_N
                  + (coor_O2 + coor_O4 + coor_O6) * mass_O
                  + (coor_C2 + coor_C4 + coor_C5 + coor_C6 + coor_C8)
                  * mass_C )
            * (1/(n_nb*mass_N + n_ob*mass_O + n_cb*mass_C));
        _P.set_coords(com_P);
        if (_S.atom_name() != "S  ")
            _S.set_coords(com_S);
        if (_B.atom_name() != "B  ")
            _B.set_coords(com_B);
    }

    /*                _     _       _
    //  _ __ ___  ___(_) __| |  ___(_)_______
    // | '__/ _ \/ __| |/ _` | / __| |_  / _ \
    // | | |  __/\__ \ | (_| | \__ \ |/ /  __/
    // |_|  \___||___/_|\__,_| |___/_/___\___|
    */
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
        _P.reset();
        _S.reset();
        _B.reset();
        _chain_type = none;
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

        _chain_type = none;
        _P.reset();
        _S.reset();
        _B.reset();
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

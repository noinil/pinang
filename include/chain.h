// -*-c++-*-

#ifndef PINANG_CHAIN_H_
#define PINANG_CHAIN_H_

#include <iostream>
#include "residue.h"

namespace pinang {
    class Chain
    {
    public:
        Chain();
        virtual ~Chain() {_residues.clear();}

        inline void reset();

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline chain_t chain_type() const;
        inline void set_chain_type(chain_t ct);

        inline Residue& m_residue(unsigned int n);
        inline int add_residue(const Residue& r);

        inline int m_chain_length() const;

        inline void pr_seq(int n) const;
        void output_cg_pos(std::ostream& o, int& n);
        void output_top_mass(std::ostream& o, int& n);
        void output_top_bond(std::ostream& o, int& n);
        void output_top_angle(std::ostream& o, int& n);
        void output_top_dihedral(std::ostream& o, int& n);
        void output_top_native(std::ostream& o);
        int m_native_contact_number();

        Chain operator+(Chain& other);

    protected:
        char _chain_ID;
        chain_t _chain_type;
        int _n_residue;
        std::vector<Residue> _residues;
    };

    /*      _           _         ___ ____
    //  ___| |__   __ _(_)_ __   |_ _|  _ \
    // / __| '_ \ / _` | | '_ \   | || | | |
    //| (__| | | | (_| | | | | |  | || |_| |
    // \___|_| |_|\__,_|_|_| |_| |___|____/
    */
    inline char Chain::chain_ID() const
    {
        return _chain_ID;
    }
    inline void Chain::set_chain_ID(char a)
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
    inline chain_t Chain::chain_type() const
    {
        return _chain_type;
    }
    inline void Chain::set_chain_type(chain_t a)
    {
        _chain_type = a;
    }

    /*
    //                _     _
    //  _ __ ___  ___(_) __| |_   _  ___
    // | '__/ _ \/ __| |/ _` | | | |/ _ \
    // | | |  __/\__ \ | (_| | |_| |  __/
    // |_|  \___||___/_|\__,_|\__,_|\___|
    */
    inline Residue& Chain::m_residue(unsigned int n)
    {
        if (_residues.empty())
        {
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~              PINANG :: Chain               ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << "ERROR: No Residues found in Chain: "
                      << _chain_ID << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n >= _residues.size())
            {
                std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
                std::cout << " ~              PINANG :: Chain               ~ " << std::endl;
                std::cout << " ============================================== " << std::endl;
                std::cerr << "ERROR: Residue index out of range in Chain: "
                          << _chain_ID << std::endl;
                exit(EXIT_SUCCESS);
            } else {
                return _residues[n];
            }
        }
    }
    inline int Chain::add_residue(const Residue& r)
    {
        // if (r.chain_ID() != _chain_ID)
        // {
        //     return 1;
        // }
        _residues.push_back(r);
        _n_residue++;
        _chain_type = r.chain_type();
        return 0;
    }

    inline int Chain::m_chain_length() const
    {
        return _n_residue;
    }

    /*             _       _
    //  _ __  _ __(_)_ __ | |_   ___  ___  __ _
    // | '_ \| '__| | '_ \| __| / __|/ _ \/ _` |
    // | |_) | |  | | | | | |_  \__ \  __/ (_| |
    // | .__/|_|  |_|_| |_|\__| |___/\___|\__, |
    // |_|                                   |_|
    */
    inline void Chain::pr_seq(int n) const
    {
        if (_chain_type == water)
        {
            return;
        }

        int i = 0;
        int j = 0;

        std::cout << " - Chain " << _chain_ID
                  << " : " << _n_residue
                  << std::endl;

        if (n == 1)
        {
            std::cout << " ";
            for (i = 0; i < _n_residue; i++) {
                // if (_residues[i].resid_name() == "HOH")
                //     continue;
                std::cout << std::setw(1) << _residues[i].short_name();
                j++;
                if (j%10 == 5)
                {
                    std::cout << " ";
                }
                if (j%10 == 0)
                {
                    std::cout << "  " << std::setw(4) << j << std::endl;
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        } else if (n == 3)
        {
            for (i = 0; i < _n_residue; i++) {
                // if (_residues[i].resid_name() == "HOH")
                //     continue;
                std::cout << std::setw(4) << _residues[i].resid_name() ;
                j++;
                if (j%10 == 0)
                {
                    std::cout << "  " << std::setw(4) << j << std::endl;
                }
            }
            std::cout << std::endl;
        }
    }

    void Chain::output_cg_pos(std::ostream& o, int& n)
    {
        if (_chain_type == water || _chain_type == other || _chain_type == none)
        {
            return;
        }
        o << " - Chain " << _chain_ID;
        o << " : " ;
        switch (_chain_type) {
        case protein:
            o << "protein";
            break;
        case DNA:
            o << "DNA";
            break;
        case RNA:
            o << "RNA";
            break;
        case water:
            o << "water";
            break;
        case ion:
            o << "ion";
            break;
        case other:
            o << "other";
            break;
        case na:
            o << "na";
            break;
        default:
            o << "unknown";
        }
        o << " : " << _n_residue
          << std::endl;

        int i = 0;
        if (_chain_type != DNA && _chain_type != RNA && _chain_type != na)
        {
            for (i = 0; i < _n_residue; i++) {
                o << std::setw(6) << ++n
                  << std::setw(5) << _residues[i].resid_name()
                  << std::setw(5) << _residues[i].resid_index() << "   "
                  << _residues[i].m_C_alpha().coordinates().x() << " "
                  << _residues[i].m_C_alpha().coordinates().y() << " "
                  << _residues[i].m_C_alpha().coordinates().z() << " "
                  << std::endl;
            }
        } else {
            _residues[0].set_cg_na();
            o << std::setw(6) << ++n
              << std::setw(5) << "S"
              << std::setw(5) << _residues[i].resid_index() << "   "
              << _residues[0].m_S().coordinates().x() << " "
              << _residues[0].m_S().coordinates().y() << " "
              << _residues[0].m_S().coordinates().z() << " "
              << std::endl;
            o << std::setw(6) << ++n
              << std::setw(5) << _residues[i].short_name()
              << std::setw(5) << _residues[i].resid_index() << "   "
              << _residues[0].m_B().coordinates().x() << " "
              << _residues[0].m_B().coordinates().y() << " "
              << _residues[0].m_B().coordinates().z() << " "
              << std::endl;
            for (i = 1; i < _n_residue; i++) {
                _residues[i].set_cg_na();
                o << std::setw(6) << ++n
                  << std::setw(5) << "P"
                  << std::setw(5) << _residues[i].resid_index() << "   "
                  << _residues[i].m_P().coordinates().x() << " "
                  << _residues[i].m_P().coordinates().y() << " "
                  << _residues[i].m_P().coordinates().z() << " "
                  << std::endl;
                o << std::setw(6) << ++n
                  << std::setw(5) << "S"
                  << std::setw(5) << _residues[i].resid_index() << "   "
                  << _residues[i].m_S().coordinates().x() << " "
                  << _residues[i].m_S().coordinates().y() << " "
                  << _residues[i].m_S().coordinates().z() << " "
                  << std::endl;
                o << std::setw(6) << ++n
                  << std::setw(5) << _residues[i].short_name()
                  << std::setw(5) << _residues[i].resid_index() << "   "
                  << _residues[i].m_B().coordinates().x() << " "
                  << _residues[i].m_B().coordinates().y() << " "
                  << _residues[i].m_B().coordinates().z() << " "
                  << std::endl;
            }
        }
        // for (i = 0; i < _n_residue; i++) {
        //     o << std::setw(6) << i+1+n
        //       << std::setw(5) << _residues[i].resid_name()
        //       << _residues[i].m_C_alpha().coordinates()
        //       << std::endl;
        // }
        o << std::endl;
    }

    void Chain::output_top_mass(std::ostream& o, int& n)
    {
        if (_chain_type == water || _chain_type == other || _chain_type == none)
        {
            return;
        }

        int i = 0;
        if (_chain_type != DNA && _chain_type != RNA && _chain_type != na)
        {
            for (i = 0; i < _n_residue; i++) {
                o << std::setw(11) << ++n
                  << std::setw(8) << _residues[i].resid_index()
                  << std::setw(8) << _residues[i].resid_name()
                  << std::setw(8) << "CA"
                  << std::setw(10) << _residues[i].resid_mass()
                  << std::setw(8)
                  << _residues[i].resid_charge()
                  << std::endl;
            }
        } else {
            o << std::setw(11) << ++n
              << std::setw(8) << _residues[0].resid_index()
              << std::setw(8) << _residues[0].resid_name()
              << std::setw(8) << "S"
              << std::setw(10) << 99.11
              << std::setw(8) << 0.0
              << std::endl;
            o << std::setw(11) << ++n
              << std::setw(8) << _residues[0].resid_index()
              << std::setw(8) << _residues[0].resid_name()
              << std::setw(8) << "B"
              << std::setw(10) << _residues[0].resid_mass()
              << std::setw(8) << 0.0
              << std::endl;
            for (i = 1; i < _n_residue; i++) {
                o << std::setw(11) << ++n
                  << std::setw(8) << _residues[i].resid_index()
                  << std::setw(8) << _residues[i].resid_name()
                  << std::setw(8) << "P"
                  << std::setw(10) << 94.93
                  << std::setw(8) << -0.6
                  << std::endl;
                o << std::setw(11) << ++n
                  << std::setw(8) << _residues[i].resid_index()
                  << std::setw(8) << _residues[i].resid_name()
                  << std::setw(8) << "S"
                  << std::setw(10) << 99.11
                  << std::setw(8) << 0.0
                  << std::endl;
                o << std::setw(11) << ++n
                  << std::setw(8) << _residues[i].resid_index()
                  << std::setw(8) << _residues[i].resid_name()
                  << std::setw(8) << "B"
                  << std::setw(10) << _residues[i].resid_mass()
                  << std::setw(8) << 0.0
                  << std::endl;
            }
        }
    }

    void Chain::output_top_bond(std::ostream& o, int& n)
    {
        if (_chain_type == water || _chain_type == other || _chain_type == none)
        {
            return;
        }

        int i = 0;
        double d = 0;
        double d_ps = 0;
        double d_sb = 0;
        double d_sp = 0;

        if (_chain_type != DNA && _chain_type != RNA && _chain_type != na)
        {
            for (i = 0; i < _n_residue - 1; i++) {
                d = resid_ca_distance(_residues[i], _residues[i+1]);
		n++;
                o << std::setw(8) << n
                  << std::setw(6) << n + 1
                  << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4)
                  << std::setw(10) << d
                  << std::setprecision(1)
                  << std::setw(8) << p_K_bond
                  << std::endl;
            }
            n++;
        } else {
            d_sb = atom_distance(_residues[0].m_S(), _residues[0].m_B());
            o << std::setw(8) << n+1
              << std::setw(6) << n+2
              << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
              << std::setw(10) << d_sb
              << std::setprecision(1) << std::setw(8) << p_K_bond
              << std::endl;
            d_sp = atom_distance(_residues[1].m_P(), _residues[0].m_S());
            o << std::setw(8) << n+1
              << std::setw(6) << n+3
              << std::setiosflags(std::ios_base::fixed)
              << std::setprecision(4)
              << std::setw(10) << d_sp
              << std::setprecision(1)
              << std::setw(8) << p_K_bond
              << std::endl;
            for (i = 1; i < _n_residue - 1; i++) {
                n += 3;
                d_ps = atom_distance(_residues[i].m_P(), _residues[i].m_S());
                o << std::setw(8) << n
                  << std::setw(6) << n+1 << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4) << std::setw(10) << d_ps
                  << std::setprecision(1) << std::setw(8) << p_K_bond
                  << std::endl;
                d_sb = atom_distance(_residues[i].m_S(), _residues[i].m_B());
                o << std::setw(8) << n+1
                  << std::setw(6) << n+2 << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4) << std::setw(10) << d_sb
                  << std::setprecision(1) << std::setw(8) << p_K_bond
                  << std::endl;
                d_sp = atom_distance(_residues[i].m_S(), _residues[i+1].m_P());
                o << std::setw(8) << n+1
                  << std::setw(6) << n+3 << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4) << std::setw(10) << d_sp
                  << std::setprecision(1) << std::setw(8) << p_K_bond
                  << std::endl;
            }
            i = _n_residue - 1;
            n += 3;
            d_ps = atom_distance(_residues[i].m_P(), _residues[i].m_S());
            o << std::setw(8) << n
              << std::setw(6) << n+1 << std::setiosflags(std::ios_base::fixed)
              << std::setprecision(4) << std::setw(10) << d_ps
              << std::setprecision(1) << std::setw(8) << p_K_bond
              << std::endl;
            d_sb = atom_distance(_residues[i].m_S(), _residues[i].m_B());
            o << std::setw(8) << n+1
              << std::setw(6) << n+2 << std::setiosflags(std::ios_base::fixed)
              << std::setprecision(4) << std::setw(10) << d_sb
              << std::setprecision(1) << std::setw(8) << p_K_bond
              << std::endl;
            n += 2;
        }
    }

    void Chain::output_top_angle(std::ostream& o, int& n)
    {
        if (_chain_type == water || _chain_type == other || _chain_type == none)
        {
            return;
        }

        int i = 0;
        double a = 0;
        Vec3d v1, v2;
        if (_chain_type != DNA && _chain_type != RNA && _chain_type != na)
        {
            for (i = 0; i < _n_residue-2; i++) {
                v1 = _residues[i].m_C_alpha().coordinates()
                    - _residues[i+1].m_C_alpha().coordinates();
                v2 = _residues[i+2].m_C_alpha().coordinates()
                    - _residues[i+1].m_C_alpha().coordinates();
                a = vec_angle_deg (v1, v2);
                o << std::setw(8) << i+1+n
                  << std::setw(6) << i+2+n
                  << std::setw(6) << i+3+n
                  << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4)
                  << std::setw(12) << a
                  << std::setprecision(1)
                  << std::setw(8) << p_K_angle
                  << std::endl;
            }
            n += _n_residue;
        } else {
            // ---------- angle BSP ----------
            v1 = _residues[0].m_B().coordinates()
                - _residues[0].m_S().coordinates();
            v2 = _residues[1].m_P().coordinates()
                - _residues[0].m_S().coordinates();
            a = vec_angle_deg (v1, v2);
            o << std::setw(8) << n+2
              << std::setw(6) << n+1
              << std::setw(6) << n+3
              << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
              << std::setw(12) << a << std::setprecision(1)
              << std::setw(8) << p_K_angle << std::endl;
            // ---------- angle SPS ----------
            v1 = _residues[0].m_S().coordinates()
                - _residues[1].m_P().coordinates();
            v2 = _residues[1].m_S().coordinates()
                - _residues[1].m_P().coordinates();
            a = vec_angle_deg (v1, v2);
            o << std::setw(8) << n+1
              << std::setw(6) << n+3
              << std::setw(6) << n+4
              << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
              << std::setw(12) << a << std::setprecision(1)
              << std::setw(8) << p_K_angle << std::endl;

            // -------------------- loop --------------------
            for (i = 1; i < _n_residue-1; i++) {
                n += 3;
                // ---------- angle PSB ----------
                v1 = _residues[i].m_P().coordinates()
                    - _residues[i].m_S().coordinates();
                v2 = _residues[i].m_B().coordinates()
                    - _residues[i].m_S().coordinates();
                a = vec_angle_deg (v1, v2);
                o << std::setw(8) << n
                  << std::setw(6) << n+1
                  << std::setw(6) << n+2
                  << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
                  << std::setw(12) << a << std::setprecision(1)
                  << std::setw(8) << p_K_angle << std::endl;
                // ---------- angle PSP ----------
                v1 = _residues[i].m_P().coordinates()
                    - _residues[i].m_S().coordinates();
                v2 = _residues[i+1].m_P().coordinates()
                    - _residues[i].m_S().coordinates();
                a = vec_angle_deg (v1, v2);
                o << std::setw(8) << n
                  << std::setw(6) << n+1
                  << std::setw(6) << n+3
                  << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
                  << std::setw(12) << a << std::setprecision(1)
                  << std::setw(8) << p_K_angle << std::endl;
                // ---------- angle BSP ----------
                v1 = _residues[i].m_B().coordinates()
                    - _residues[i].m_S().coordinates();
                v2 = _residues[i+1].m_P().coordinates()
                    - _residues[i].m_S().coordinates();
                a = vec_angle_deg (v1, v2);
                o << std::setw(8) << n+2
                  << std::setw(6) << n+1
                  << std::setw(6) << n+3
                  << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
                  << std::setw(12) << a << std::setprecision(1)
                  << std::setw(8) << p_K_angle << std::endl;
                // ---------- angle SPS ----------
                v1 = _residues[i].m_S().coordinates()
                    - _residues[i+1].m_P().coordinates();
                v2 = _residues[i+1].m_S().coordinates()
                    - _residues[i+1].m_P().coordinates();
                a = vec_angle_deg (v1, v2);
                o << std::setw(8) << n+1
                  << std::setw(6) << n+3
                  << std::setw(6) << n+4
                  << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
                  << std::setw(12) << a << std::setprecision(1)
                  << std::setw(8) << p_K_angle << std::endl;
            }
            n += 3;
            i = _n_residue - 1;
            // ---------- angle PSB ----------
            v1 = _residues[i].m_P().coordinates()
                - _residues[i].m_S().coordinates();
            v2 = _residues[i].m_B().coordinates()
                - _residues[i].m_S().coordinates();
            a = vec_angle_deg (v1, v2);
            o << std::setw(8) << n
              << std::setw(6) << n+1
              << std::setw(6) << n+2
              << std::setiosflags(std::ios_base::fixed) << std::setprecision(4)
              << std::setw(12) << a << std::setprecision(1)
              << std::setw(8) << p_K_angle << std::endl;
            n += 2;
        }
    }

    void Chain::output_top_dihedral(std::ostream& o, int& n)
    {
        if (_chain_type == water || _chain_type == other || _chain_type == none)
        {
            return;
        }

        int i = 0;
        double d = 0;           // dihedral
        Vec3d v1, v2, v3, n1, n2;
        if (_chain_type != DNA && _chain_type != RNA && _chain_type != na)
        {
            for (i = 0; i < _n_residue-3; i++) {
                v1 = _residues[i].m_C_alpha().coordinates()
                    - _residues[i+1].m_C_alpha().coordinates();
                v2 = _residues[i+2].m_C_alpha().coordinates()
                    - _residues[i+1].m_C_alpha().coordinates();
                v3 = _residues[i+2].m_C_alpha().coordinates()
                    - _residues[i+3].m_C_alpha().coordinates();
                n1 = v1 ^ v2;
                n2 = v2 ^ v3;
                d = vec_angle_deg (n1, n2);
                o << std::setw(8) << i+1+n
                  << std::setw(6) << i+2+n
                  << std::setw(6) << i+3+n
                  << std::setw(6) << i+4+n
                  << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4)
                  << std::setw(12) << d
                  << std::setprecision(1)
                  << std::setw(8) << p_K_dihedral_1
                  << std::setw(8) << p_K_dihedral_3
                  << std::endl;
            }
            n += _n_residue;
        } else {
            v1 = _residues[0].m_S().coordinates()
                - _residues[1].m_P().coordinates();
            v2 = _residues[1].m_S().coordinates()
                - _residues[1].m_P().coordinates();
            v3 = _residues[1].m_S().coordinates()
                - _residues[2].m_P().coordinates();
            n1 = v1 ^ v2;
            n2 = v2 ^ v3;
            d = vec_angle_deg (n1, n2);
            o << std::setw(8) << n + 1
              << std::setw(6) << n + 3
              << std::setw(6) << n + 4
              << std::setw(6) << n + 6
              << std::setiosflags(std::ios_base::fixed)
              << std::setprecision(4)
              << std::setw(12) << d
              << std::setprecision(1)
              << std::setw(8) << p_K_dihedral_1
              << std::setw(8) << p_K_dihedral_3
              << std::endl;

            for (i = 1; i < _n_residue-1; i++) {
                n += 3;
                // ---------- PSPS ----------
                v1 = _residues[i].m_P().coordinates()
                    - _residues[i].m_S().coordinates();
                v2 = _residues[i+1].m_P().coordinates()
                    - _residues[i].m_S().coordinates();
                v3 = _residues[i+1].m_P().coordinates()
                    - _residues[i+1].m_S().coordinates();
                n1 = v1 ^ v2;
                n2 = v2 ^ v3;
                d = vec_angle_deg (n1, n2);
                o << std::setw(8) << n
                  << std::setw(6) << n + 1
                  << std::setw(6) << n + 3
                  << std::setw(6) << 4 + n
                  << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4)
                  << std::setw(12) << d
                  << std::setprecision(1)
                  << std::setw(8) << p_K_dihedral_1
                  << std::setw(8) << p_K_dihedral_3
                  << std::endl;

                if (i == _n_residue - 2)
                    break;
                // ---------- SPSP ----------
                v1 = _residues[i].m_S().coordinates()
                    - _residues[i+1].m_P().coordinates();
                v2 = _residues[i+1].m_S().coordinates()
                    - _residues[i+1].m_P().coordinates();
                v3 = _residues[i+1].m_S().coordinates()
                    - _residues[i+2].m_P().coordinates();
                n1 = v1 ^ v2;
                n2 = v2 ^ v3;
                d = vec_angle_deg (n1, n2);
                o << std::setw(8) << n + 1
                  << std::setw(6) << n + 3
                  << std::setw(6) << n + 4
                  << std::setw(6) << n + 6
                  << std::setiosflags(std::ios_base::fixed)
                  << std::setprecision(4)
                  << std::setw(12) << d
                  << std::setprecision(1)
                  << std::setw(8) << p_K_dihedral_1
                  << std::setw(8) << p_K_dihedral_3
                  << std::endl;
            }
            n += 5;
        }
    }

    void Chain::output_top_native(std::ostream& o)
    {
        int i = 0, j = 0;
        double d = -1, f = -1;
        chain_t cti, ctj;
        for (i = 0; i < _n_residue-4; i++) {
            cti = _residues[i].chain_type();
            if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
                continue;
            for (j = i + 4; j < _n_residue; j++) {
                ctj = _residues[j].chain_type();
                if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
                    continue;
                d = resid_min_distance(_residues[i], _residues[j]);
                if ( d < g_cutoff)
                {
                    f = resid_ca_distance(_residues[i], _residues[j]);
                    o << std::setw(8) << i+1
                      << std::setw(6) << j+1
                      << std::setiosflags(std::ios_base::fixed)
                      << std::setprecision(2)
                      << std::setw(8) << p_K_native
                      << std::setprecision(4)
                      << std::setw(10) << f
                      << std::endl;
                }
            }
        }
    }

    int Chain::m_native_contact_number()
    {
        int i = 0, j = 0;
        double d = -1;
        int n = 0;
        chain_t cti, ctj;

        for (i = 0; i < _n_residue-4; i++) {
            cti = _residues[i].chain_type();
            if (cti == water || cti == DNA || cti == RNA || cti == na || cti == ion)
                continue;
            for (j = i + 4; j < _n_residue; j++) {
                ctj = _residues[j].chain_type();
                if (ctj == water || ctj == DNA || ctj == RNA || ctj == na || ctj == ion)
                    continue;
                d = resid_min_distance(_residues[i], _residues[j]);
                if ( d < g_cutoff)
                    n++;
            }
        }
        return n;
    }

// Chain -------------------------------------------------------------------
    inline Chain::Chain()
    {
        _chain_ID = -1;
        _chain_type = none;
        _residues.clear();
        _n_residue = 0;
    }

    inline void Chain::reset()
    {
        _chain_ID = -1;
        _chain_type = none;
        _residues.clear();
        _n_residue = 0;
    }

    Chain Chain::operator+(Chain& other)
    {
        int i = 0;
        Chain c1;
        c1.set_chain_ID(other._chain_ID);
        if (_n_residue > 0)
        {
            for (i = 0; i < _n_residue; i++) {
                c1._residues.push_back(_residues[i]);
                c1._n_residue++;
            }
        }
        for (i = 0; i < other._n_residue; i++) {
            c1._residues.push_back(other._residues[i]);
            c1._n_residue++;
        }

        return c1;
    }

    inline std::ostream& operator<<(std::ostream& o, Chain& c)
    {
        // o << "Chain "
        //   << std::setw(4) << c.chain_ID() << ":  "
        //   << std::endl;
        int i = 0;
        for (i = 0; i < c.m_chain_length(); i++) {
            o << c.m_residue(i) ;
        }
        o << "TER   " << std::endl;
        return o;
    }

}

#endif

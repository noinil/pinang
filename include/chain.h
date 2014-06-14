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

        inline Residue& m_residue(unsigned int n);
        inline int add_residue(const Residue& r);

        inline int m_chain_length() const;

        inline void pr_seq(int n) const;
        void output_ca_pos(std::ostream& o, int n);
        void output_top_mass(std::ostream& o, int n);
        void output_top_bond(std::ostream& o, int n);
        void output_top_angle(std::ostream& o, int n);
        void output_top_dihedral(std::ostream& o, int n);
        void output_top_native(std::ostream& o);

        Chain operator+(Chain& other);

    protected:
        char _chain_ID;
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
            if (n < 0 || n >= _residues.size())
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
        if (_residues[0].resid_name() == "HOH")
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

    void Chain::output_ca_pos(std::ostream& o, int n)
    {
        if (_residues[0].resid_name() == "HOH")
        {
            return;
        }
        o << " - Chain " << _chain_ID
          << " : " << _n_residue
          << std::endl;

        int i = 0;
        for (i = 0; i < _n_residue; i++) {
            o << std::setw(6) << i+1+n
              << std::setw(5) << _residues[i].resid_name()
              << _residues[i].m_C_alpha().coordinates()
              << std::endl;
        }
        o << std::endl;
    }

    void Chain::output_top_mass(std::ostream& o, int n)
    {
        if (_residues[0].resid_name() == "HOH")
        {
            return;
        }

        int i = 0;
        for (i = 0; i < _n_residue; i++) {
                o << std::setw(11) << i+1+n
                  << std::setw(10) << _residues[i].resid_mass()
                  << std::setw(8)
                  << _residues[i].resid_charge()
                  << std::endl;
        }
    }

    void Chain::output_top_bond(std::ostream& o, int n)
    {
        if (_residues[0].resid_name() == "HOH")
        {
            return;
        }

        int i = 0;
        for (i = 0; i < _n_residue - 1; i++) {
            o << std::setw(8) << i+1+n
              << std::setw(6) << i+2+n
              << std::setw(8) << p_K_bond
              << std::endl;
        }
    }

    void Chain::output_top_angle(std::ostream& o, int n)
    {
        if (_residues[0].resid_name() == "HOH")
        {
            return;
        }

        int i = 0;
        for (i = 0; i < _n_residue-2; i++) {
            o << std::setw(8) << i+1+n
              << std::setw(6) << i+2+n
              << std::setw(6) << i+3+n
              << std::setw(8) << p_K_angle
              << std::endl;
        }
    }

    void Chain::output_top_dihedral(std::ostream& o, int n)
    {
        if (_residues[0].resid_name() == "HOH")
        {
            return;
        }

        int i = 0;
        for (i = 0; i < _n_residue-3; i++) {
            o << std::setw(8) << i+1+n
              << std::setw(6) << i+2+n
              << std::setw(6) << i+3+n
              << std::setw(6) << i+4+n
              << std::setw(8) << p_K_dihedral_1
              << std::setw(8) << p_K_dihedral_3
              << std::endl;
        }
    }

    void Chain::output_top_native(std::ostream& o)
    {
        int i = 0, j = 0;
        double d = -1, f = -1;
        for (i = 0; i < _n_residue-4; i++) {
            if (_residues[i].resid_name() == "HOH")
                continue;
            for (j = i + 4; j < _n_residue; j++) {
                if (_residues[j].resid_name() == "HOH")
                    continue;
                d = resid_min_distance(_residues[i], _residues[j]);
                if ( d < g_cutoff)
                {
                    f = resid_ca_distance(_residues[i], _residues[j]);
                    o << std::setw(8) << i+1
                      << std::setw(6) << j+1
                      << std::setw(8) << p_K_native
                      << std::setiosflags(std::ios_base::fixed)
                      << std::setprecision(2)
                      << std::setw(8) << f
                      << std::setw(8) << d
                      << std::endl;
                }
            }
        }
    }

// Chain -------------------------------------------------------------------
    inline Chain::Chain()
    {
        _chain_ID = -1;
        _residues.clear();
        _n_residue = 0;
    }

    inline void Chain::reset()
    {
        _chain_ID = -1;
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

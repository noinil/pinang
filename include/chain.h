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
        void output_ca_pos(std::ostream& o);

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
        int i = 0;
        int j = 0;

        if (n == 1)
        {
            std::cout << " ";
            for (i = 0; i < _n_residue; i++) {
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

    void Chain::output_ca_pos(std::ostream& o)
    {
        int i = 0;
        for (i = 0; i < _n_residue; i++) {
            if (_residues[i].resid_name() != "HOH")
            {
                o << std::setw(6) << _residues[i].resid_index()
                  << std::setw(5) << _residues[i].resid_name()
                  << _residues[i].m_C_alpha().coordinates()
                  << std::endl;
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

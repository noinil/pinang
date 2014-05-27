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

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline Residue& m_residue(unsigned int n);
        inline void add_residue(const Residue& r);

        inline int m_chain_length() const;

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
            std::cerr << "ERROR: No Residues found in Chain: "
                      << _chain_ID << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _residues.size())
            {
                std::cerr << "ERROR: Residue index out of range in Chain: "
                          << _chain_ID << std::endl;
                exit(EXIT_SUCCESS);
            } else {
                return _residues[n];
            }
        }
    }
    inline void Chain::add_residue(const Residue& r)
    {
        _residues.push_back(r);
        _n_residue++;
    }

    inline int m_chain_length() const
    {
        return _n_residue;
    }

    // Chain -------------------------------------------------------------------
    inline Chain::Chain()
    {
        _chain_ID = 0;
        _residues.clear();
        _n_residue = 0;
    }

}

#endif

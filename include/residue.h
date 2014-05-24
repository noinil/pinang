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

        inline std::string resid_name() const;
        inline void set_resid_name(const std::string& s);

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline unsigned int resid_index() const;
        inline void set_resid_index(unsigned int i);

        inline const Atom& m_atom(unsigned int n) const;
        inline void add_atom(const Atom& a);

        inline const Vec3d& pos_Ca() const;
        inline void set_pos_Ca();
        inline void update_pos_Ca(double dt);

        inline const Vec3d& vel_Ca() const;
        inline void update_vel_Ca(double dt, const Vec3d& acceleration);

    protected:
        std::string _resid_name;
        char _chain_ID;
        unsigned int _resid_index;
        std::vector<Atom> _atoms;
        Vec3d _pos_Ca;          // position of C_alpha;
        Vec3d _vel_Ca;          // velocity of C_alpha;
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
    inline const Atom& Residue::m_atom(unsigned int n) const
    {
        if (_atoms.empty())
        {
            std::cout << "No Atoms found in Residue: " << _resid_index << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _atoms.size())
            {
                std::cout << "Atom index out of range in Residue: "
                          << _resid_index << std::endl;
                exit(EXIT_SUCCESS);
            } else {
                return _atoms[n];
            }
        }
    }
    inline void Residue::add_atom(const Atom& a)
    {
        _atoms.push_back(a);
    }

    // Residue------------------------------------------------------------------
    inline Residue::Residue()
    {
        _resid_name = "";
        _chain_ID = 0;
        _resid_index = 0;
        _atoms.clear();
        _pos_Ca = Vec3d(0, 0, 0);
        _vel_Ca = Vec3d(0, 0, 0);
    }

}

#endif

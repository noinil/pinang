// -*-c++-*-

#ifndef BINANG_RESIDUE_H
#define BINANG_RESIDUE_H

#include <iostream>
#include "atom.h"

namespace binang {
    class Residue
    {
    public:
        Residue();
        virtual ~Residue();

        inline std::string resid_name() const;
        inline void set_resid_name(const std::string& s);

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline unsigned int resid_index() const;
        inline void set_resid_index(unsigned int i);

        inline const Atom& m_atom(int n) const;
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


}

#endif

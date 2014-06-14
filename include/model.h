// -*-c++-*-

#ifndef PINANG_MODEL_H_
#define PINANG_MODEL_H_

#include <iostream>
#include "chain.h"

namespace pinang {

    class Model
    {
    public:
        Model();
        virtual ~Model() {_chains.clear();}

        inline void reset();

        inline int model_ID() const;
        inline void set_model_ID(int n);

        inline Chain& m_chain(unsigned int n);
        inline void add_chain(const Chain& c);

        inline void print_sequence(int n) const;

        inline int m_model_size() const;

        void output_ca_pos(std::ostream& o);
        void output_top_mass(std::ostream& o);
        void output_top_bond(std::ostream& o);
        void output_top_angle(std::ostream& o);
        void output_top_dihedral(std::ostream& o);

        void output_top_native(std::ostream& o);

    protected:
        int _model_ID;
        std::vector<Chain> _chains;
        int _n_chain;
    };

    inline int Model::model_ID() const
    {
        return _model_ID;
    }
    inline void Model::set_model_ID(int n)
    {
        _model_ID = n;
    }

    /*       _           _
    //   ___| |__   __ _(_)_ __
    //  / __| '_ \ / _` | | '_ \
    // | (__| | | | (_| | | | | |
    //  \___|_| |_|\__,_|_|_| |_|
    */
    inline Chain& Model::m_chain(unsigned int n)
    {
        if (_chains.empty())
        {
            std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
            std::cout << " ~              PINANG :: Model               ~ " << std::endl;
            std::cout << " ============================================== " << std::endl;
            std::cerr << "ERROR: No Chains found in Model: "
                      << _model_ID << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _chains.size())
            {
                std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
                std::cout << " ~              PINANG :: Model               ~ " << std::endl;
                std::cout << " ============================================== " << std::endl;
                std::cerr << "ERROR: Chain number out of range in Model: "
                          << _model_ID << std::endl;
                exit(EXIT_SUCCESS);
            } else {
                return _chains[n];
            }
        }
    }
    inline void Model::add_chain(const Chain& c)
    {
        _chains.push_back(c);
        _n_chain++;
    }

    /*             _       _
    //  _ __  _ __(_)_ __ | |_   ___  ___  __ _
    // | '_ \| '__| | '_ \| __| / __|/ _ \/ _` |
    // | |_) | |  | | | | | |_  \__ \  __/ (_| |
    // | .__/|_|  |_|_| |_|\__| |___/\___|\__, |
    // |_|                                   |_|
    */
    inline void Model::print_sequence(int n) const
    {
        int i = 0;
        for (i = 0; i < _n_chain; i++) {
            _chains[i].pr_seq(n);
        }

    }

    inline int Model::m_model_size() const
    {
        return _n_chain;
    }

    // Model ===================================================================
    inline Model::Model()
    {
        _model_ID = 0;
        _chains.clear();
        _n_chain = 0;
    }

    inline void Model::reset()
    {
        _model_ID = 0;
        _chains.clear();
        _n_chain = 0;
    }


    /*              _               _     _                    _
    //   ___  _   _| |_ _ __  _   _| |_  | |_ ___  _ __   ___ | | ___   __ _ _   _
    //  / _ \| | | | __| '_ \| | | | __| | __/ _ \| '_ \ / _ \| |/ _ \ / _` | | | |
    // | (_) | |_| | |_| |_) | |_| | |_  | || (_) | |_) | (_) | | (_) | (_| | |_| |
    //  \___/ \__,_|\__| .__/ \__,_|\__|  \__\___/| .__/ \___/|_|\___/ \__, |\__, |
    //                 |_|                        |_|                  |___/ |___/
    */
    void Model::output_ca_pos(std::ostream& o)
    {
        int i = 0;
        int n = 0;
        for (i = 0; i < _n_chain; i++) {
            if (_chains[i].m_residue(0).resid_name() == "HOH")
                continue;
            _chains[i].output_ca_pos(o, n);
            n += _chains[i].m_chain_length();
        }
    }

    void Model::output_top_mass(std::ostream& o)
    {
        int i = 0;
        int n = 0;
        o << "[ particles ]" << std::endl;
        o << "# "
          << std::setw(9) << "index"
          << std::setw(10) << "mass"
          << std::setw(8) << "charge"
          << std::endl;

        for (i = 0; i < _n_chain; i++) {
            if (_chains[i].m_residue(0).resid_name() == "HOH")
                continue;
            _chains[i].output_top_mass(o, n);
            n += _chains[i].m_chain_length();
        }
        o << std::endl;
    }

    void Model::output_top_bond(std::ostream& o)
    {
        int i = 0;
        int n = 0;
        o << "[ bonds ]" << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(8) << "K_b"
          << std::endl;

        for (i = 0; i < _n_chain; i++) {
            if (_chains[i].m_residue(0).resid_name() == "HOH")
                continue;
            _chains[i].output_top_bond(o, n);
            n += _chains[i].m_chain_length();
        }
        o << std::endl;
    }

    void Model::output_top_angle(std::ostream& o)
    {
        int i = 0;
        int n = 0;
        o << "[ angles ]" << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(6) << "pk"
          << std::setw(8) << "K_a"
          << std::endl;

        for (i = 0; i < _n_chain; i++) {
            if (_chains[i].m_residue(0).resid_name() == "HOH")
                continue;
            _chains[i].output_top_angle(o, n);
            n += _chains[i].m_chain_length();
        }
        o << std::endl;
    }

    void Model::output_top_dihedral(std::ostream& o)
    {
        int i = 0;
        int n = 0;
        o << "[ dihedrals ]" << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(6) << "pk"
          << std::setw(6) << "pl"
          << std::setw(8) << "K_d_1"
          << std::setw(8) << "K_d_3"
          << std::endl;

        for (i = 0; i < _n_chain; i++) {
            if (_chains[i].m_residue(0).resid_name() == "HOH")
                continue;
            _chains[i].output_top_dihedral(o, n);
            n += _chains[i].m_chain_length();
        }
        o << std::endl;
    }

    void Model::output_top_native(std::ostream& o)
    {
        int i = 0;
        Chain c0;
        o << "[ native ]" << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(8) << "eps"
          << std::setw(10) << "sigma"
          << std::endl;

        for (i = 0; i < _n_chain; i++) {
            if (_chains[i].m_residue(0).resid_name() == "HOH")
                continue;
            c0 = c0 + _chains[i];
        }

        c0.output_top_native(o);
        o << std::endl;
    }

    /*            _
    //   ___  ___| |_ _ __ ___  __ _ _ __ ___
    //  / _ \/ __| __| '__/ _ \/ _` | '_ ` _ \
    // | (_) \__ \ |_| | |  __/ (_| | | | | | |
    //  \___/|___/\__|_|  \___|\__,_|_| |_| |_|
    */
    inline std::ostream& operator<<(std::ostream& o, Model& m)
    {
        o << "MODEL "
          << std::setw(8) << m.model_ID()
          << std::endl;
        int i = 0;
        for (i = 0; i < m.m_model_size(); i++) {
            o << m.m_chain(i) ;
        }
        o << "ENDMDL" << std::endl;
        return o;
    }

}

#endif

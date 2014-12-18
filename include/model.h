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
        inline void output_fasta(std::ostream & f_fasta, std::string s) const;

        inline int m_model_size() const;

        void output_cg_pos(std::ostream& o);

        void output_top_mass(std::ostream& o);
        void output_top_bond(std::ostream& o);
        void output_top_angle(std::ostream& o);
        void output_top_dihedral(std::ostream& o);

        void output_top_nonbonded(std::ostream& o);

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
            if (n >= _chains.size())
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
    inline void Model::output_fasta(std::ostream & f_fasta, std::string s) const
    {
        int i = 0;
        for (i = 0; i < _n_chain; i++) {
            _chains[i].output_fasta(f_fasta, s);
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
    void Model::output_cg_pos(std::ostream& o)
    {
        int i = 0;
        int n = 0;
        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            _chains[i].output_cg_pos(o, n);
        }
    }

    void Model::output_top_mass(std::ostream& o)
    {
        int i = 0;
        int n = 0;

        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            else if (ct == DNA || ct == RNA || ct == na)
                n += _chains[i].m_chain_length() * 3 - 1;
            else
                n += _chains[i].m_chain_length();
        }

        o << "[ particles ]"
          << std::setw(6) << n
          << std::endl;
        o << "# "
          << std::setw(9) << "index"
          << std::setw(8) << "resid"
          << std::setw(8) << "resname"
          << std::setw(8) << "atom"
          << std::setw(10) << "mass"
          << std::setw(8) << "charge"
          << std::endl;

        n = 0;
        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            _chains[i].output_top_mass(o, n);
        }
        o << std::endl;
    }

    void Model::output_top_bond(std::ostream& o)
    {
        int i = 0;
        int n = 0;

        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            else if (ct == DNA || ct == RNA || ct == na)
                n += _chains[i].m_chain_length() * 3 - 2;
            else
                n += _chains[i].m_chain_length() - 1;
        }

        o << "[ bonds ]"
          << std::setw(6) << n
          << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(10) << "r0"
          << std::setw(8) << "K_b"
          << std::endl;

        n = 0;
        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            _chains[i].output_top_bond(o, n);
        }
        o << std::endl;
    }

    void Model::output_top_angle(std::ostream& o)
    {
        int i = 0;
        int n = 0;

        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            else if (ct == DNA || ct == RNA || ct == na)
                n += _chains[i].m_chain_length() * 4 - 5;
            else
                n += _chains[i].m_chain_length()-2;
        }

        o << "[ angles ]"
          << std::setw(6) << n
          << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(6) << "pk"
          << std::setw(12) << "theta_0"
          << std::setw(8) << "K_a"
          << std::endl;

        n = 0;
        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            _chains[i].output_top_angle(o, n);
        }
        o << std::endl;
    }

    void Model::output_top_dihedral(std::ostream& o)
    {
        int i = 0;
        int n = 0;

        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            else if (ct == DNA || ct == RNA || ct == na)
                n += _chains[i].m_chain_length() * 2 - 4;
            else
                n += _chains[i].m_chain_length()-3;
        }

        o << "[ dihedrals ]"
          << std::setw(6) << n
          << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(6) << "pk"
          << std::setw(6) << "pl"
          << std::setw(12) << "phi_0"
          << std::setw(8) << "K_d_1"
          << std::setw(8) << "K_d_3"
          << std::endl;

        n = 0;
        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            _chains[i].output_top_dihedral(o, n);
        }
        o << std::endl;
    }

    void Model::output_top_nonbonded(std::ostream& o)
    {
        int i = 0;
        Chain c0;
        Chain c_tmp;
        Residue r_tmp;

        for (i = 0; i < _n_chain; i++) {
            chain_t ct = _chains[i].chain_type();
            if (ct == water || ct == other || ct == none)
                continue;
            if (ct == DNA || ct == RNA || ct == na)
            {
                c_tmp.reset();

                r_tmp.reset();
                r_tmp.add_atom(_chains[i].m_residue(0).m_S());
                r_tmp.set_chain_type(ct);
                c_tmp.add_residue(r_tmp);

                r_tmp.reset();
                r_tmp.add_atom(_chains[i].m_residue(0).m_B());
                r_tmp.set_chain_type(ct);
                c_tmp.add_residue(r_tmp);

                for (int j = 1; j < _chains[i].m_chain_length(); j++) {
                    r_tmp.reset();
                    r_tmp.add_atom(_chains[i].m_residue(j).m_P());
                    r_tmp.set_chain_type(ct);
                    c_tmp.add_residue(r_tmp);
                    r_tmp.reset();
                    r_tmp.add_atom(_chains[i].m_residue(j).m_S());
                    r_tmp.set_chain_type(ct);
                    c_tmp.add_residue(r_tmp);
                    r_tmp.reset();
                    r_tmp.add_atom(_chains[i].m_residue(j).m_B());
                    r_tmp.set_chain_type(ct);
                    c_tmp.add_residue(r_tmp);
                }
                c_tmp.set_chain_type(ct);
                c0 = c0 + c_tmp;
                continue;
            }
            c0 = c0 + _chains[i];
        }

        o << "[ native ]"
          << std::setw(6) << c0.m_native_contact_number()
          << std::endl;
        o << "# "
          << std::setw(6) << "pi"
          << std::setw(6) << "pj"
          << std::setw(8) << "eps"
          << std::setw(10) << "sigma"
          << std::endl;

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

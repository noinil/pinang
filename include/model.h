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
            std::cerr << "ERROR: No Chains found in Model: "
                      << _model_ID << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            if (n < 0 || n >= _chains.size())
            {
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
            std::cout << " - Chain " << _chains[i].chain_ID() << std::endl;
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

// -*-c++-*-

#ifndef PINANG_PARTICLE_H_
#define PINANG_PARTICLE_H_

#include "vec3d.h"

#include <string>
#include <sstream>

namespace pinang {
    class Particle
    {
    public:
        inline Particle();
        virtual ~Particle() {};

        inline void reset();

        inline std::string atom_name() const;
        inline void set_atom_name(const std::string& s);

        inline std::string resid_name() const;
        inline void set_resid_name(const std::string& s);

        inline int resid_index() const;
        inline void set_resid_index(int i);

        inline double charge() const;
        inline void set_charge(double c);

        inline double mass() const;
        inline void set_mass(double m);

    protected:
        std::string _atom_name;
        std::string _resid_name;
        int _resid_index;
        double _charge;
        double _mass;
    };

    /*        _
    //   __ _| |_ ___  _ __ ___    _ __   __ _ _ __ ___   ___
    //  / _` | __/ _ \| '_ ` _ \  | '_ \ / _` | '_ ` _ \ / _ \
    // | (_| | || (_) | | | | | | | | | | (_| | | | | | |  __/
    //  \__,_|\__\___/|_| |_| |_| |_| |_|\__,_|_| |_| |_|\___|
    */
    inline std::string Particle::atom_name() const
    {
        return _atom_name;
    }
    inline void Particle::set_atom_name(const std::string& s)
    {
        size_t sz = 3;
        _atom_name = s;

        if (_atom_name.size() < sz)
        {
            _atom_name.resize(sz, ' ');
        }
    }

    /*                _     _
    //  _ __ ___  ___(_) __| |  _ __   __ _ _ __ ___   ___
    // | '__/ _ \/ __| |/ _` | | '_ \ / _` | '_ ` _ \ / _ \
    // | | |  __/\__ \ | (_| | | | | | (_| | | | | | |  __/
    // |_|  \___||___/_|\__,_| |_| |_|\__,_|_| |_| |_|\___|
    */
    inline std::string Particle::resid_name() const
    {
        return _resid_name;
    }
    inline void Particle::set_resid_name(const std::string& s)
    {
        _resid_name = s;
    }

    /*                _     _   _           _
    //  _ __ ___  ___(_) __| | (_)_ __   __| | _____  __
    // | '__/ _ \/ __| |/ _` | | | '_ \ / _` |/ _ \ \/ /
    // | | |  __/\__ \ | (_| | | | | | | (_| |  __/>  <
    // |_|  \___||___/_|\__,_| |_|_| |_|\__,_|\___/_/\_\
    */
    inline int Particle::resid_index() const
    {
        return _resid_index;
    }
    inline void Particle::set_resid_index(int i)
    {
        _resid_index = i;
    }

    /*       _
    //   ___| |__   __ _ _ __ __ _  ___
    //  / __| '_ \ / _` | '__/ _` |/ _ \
    // | (__| | | | (_| | | | (_| |  __/
    //  \___|_| |_|\__,_|_|  \__, |\___|
    //                       |___/
    */
    inline double Particle::charge() const
    {
        return _charge;
    }
    inline void Particle::set_charge(double c)
    {
        _charge = c;
    }

    inline double Particle::mass() const
    {
        return _mass;
    }
    inline void Particle::set_mass(double m)
    {
        _mass = m;
    }

    // Particle::Particle ==============================================================
    inline Particle::Particle()
    {
        _atom_name = "";
        _resid_name = "";
        _resid_index = 0;
        _charge = 0;
        _mass = 0;
    }

    inline void Particle::reset()
    {
        _atom_name = "";
        _resid_name = "";
        _resid_index = 0;
        _charge = 0;
        _mass = 0;
    }

    /*   ___        _              _____                 _   _
    //  / _ \ _   _| |_ ___ _ __  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
    // | | | | | | | __/ _ \ '__| | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
    // | |_| | |_| | ||  __/ |    |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
    //  \___/ \__,_|\__\___|_|    |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
    */
    std::istream& operator>>(std::istream& i, Particle& p)
    {
        std::istringstream tmp_sstr;
        std::string top_line;

        std::string tmp_s;
        int tmp_i;
        double tmp_d;

        std::getline(i, top_line);
        tmp_sstr.str (top_line);

        tmp_sstr >> tmp_i;
        std::cout << ".";

        tmp_sstr >> tmp_i;
        p.set_resid_index(tmp_i);

        tmp_sstr >> tmp_s;
        p.set_resid_name(tmp_s);

        tmp_sstr >> tmp_s;
        p.set_atom_name(tmp_s);

        tmp_sstr >> tmp_d;
        p.set_mass(tmp_d);

        tmp_sstr >> tmp_d;
        p.set_charge(tmp_d);

        if (!i) return i;
        return i;
    }

}
#endif

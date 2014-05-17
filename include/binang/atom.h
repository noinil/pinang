// -*-c++-*-

#include <iostream>
#include "vec3d.h"

namespace binang {
    class Atom
    {
    public:
        inline Atom();
        // virtual ~Atom();

        inline unsigned int serial() const;
        inline void set_serial(unsigned int i);

        inline double occupancy() const ;
        inline void set_occupancy(double o);
        inline double temperature_factor() const;
        inline void set_temperature_factor(double f);

        inline std::string segment_id() const;
        inline void set_segment_id(const std::string s);

        inline const Vec3d &coordinates() const;
        inline void set_coords(const Vec3d &coors);
    protected:
        unsigned int _serial;
        std::string _atom_name;
        char _alt_loc;
        std::string _resid_name;
        char _chain_ID;
        unsigned int _res_index;
        char _insert_code;
        Vec3d _coordinate;
        Vec3d _velocity;
        // Vec3d _accelaration;
        double _occupancy;
        double _temp_factor;
        std::string _seg_ID;
        std::string _element;
        std::string _charge;
    };

    inline unsigned int Atom::serial() const
    {
        return _serial;
    }
    inline void Atom::set_serial(unsigned int i)
    {
        _serial = i;
    }

    inline double Atom::occupancy() const
    {
        return _occupancy;
    }
    inline void Atom::set_occupancy(double o)
    {
        _occupancy = o;
    }
    inline double Atom::temperature_factor() const
    {
        return _temp_factor;
    }
    inline void Atom::set_temperature_factor(double f)
    {
        _temp_factor = f;
    }

    inline std::string Atom::segment_id() const
    {
        return _seg_ID;
    }
    inline void Atom::set_segment_id(const std::string s)
    {
        _seg_ID = s;
    }

    inline const Vec3d & Atom::coordinates() const
    {
        return _coordinate;
    }
    inline void Atom::set_coords(const Vec3d &coors)
    {
        _coordinate = coors;
    }

    inline Atom::Atom()
    {
        _occupancy = 0;
    }
}

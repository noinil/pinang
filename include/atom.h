// -*-c++-*-

#ifndef PINANG_ATOM_H_
#define PINANG_ATOM_H_

#include "vec3d.h"

#include <string>
#include <sstream>

namespace pinang {
    class Atom
    {
    public:
        inline Atom();
        virtual ~Atom() {};

        inline void reset();

        inline std::string atom_flag() const;
        inline void set_atom_flag(const std::string& s);

        inline unsigned int serial() const;
        inline void set_serial(unsigned int i);

        inline std::string atom_name() const;
        inline void set_atom_name(const std::string& s);

        inline char alt_loc() const;
        inline void set_alt_loc(char a);

        inline std::string resid_name() const;
        inline void set_resid_name(const std::string& s);

        inline char chain_ID() const;
        inline void set_chain_ID(char a);

        inline unsigned int resid_index() const;
        inline void set_resid_index(unsigned int i);

        inline char icode() const;
        inline void set_icode(char a);

        inline const Vec3d& coordinates() const;
        inline void set_coords(const Vec3d& coors);

        inline const Vec3d& velocities() const;
        inline void set_velocities(const Vec3d& velos);

        inline double occupancy() const ;
        inline void set_occupancy(double o);

        inline double temperature_factor() const;
        inline void set_temperature_factor(double f);

        inline std::string segment_ID() const;
        inline void set_segment_ID(const std::string& s);

        inline std::string element() const;
        inline void set_element(const std::string& s);

        inline std::string charge() const;
        inline void set_charge(const std::string& s);

    protected:
        std::string _atom_flag;
        unsigned int _serial;
        std::string _atom_name;
        char _alt_loc;
        std::string _resid_name;
        char _chain_ID;
        unsigned int _resid_index;
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

    inline std::string Atom::atom_flag() const
    {
        return _atom_flag;
    }
    inline void Atom::set_atom_flag(const std::string& s)
    {
        _atom_flag = s;
    }

    /*                _       _
    //  ___  ___ _ __(_) __ _| |
    // / __|/ _ \ '__| |/ _` | |
    // \__ \  __/ |  | | (_| | |
    // |___/\___|_|  |_|\__,_|_|
    */
    inline unsigned int Atom::serial() const
    {
        return _serial;
    }
    inline void Atom::set_serial(unsigned int i)
    {
        _serial = i;
    }

    /*        _
    //   __ _| |_ ___  _ __ ___    _ __   __ _ _ __ ___   ___
    //  / _` | __/ _ \| '_ ` _ \  | '_ \ / _` | '_ ` _ \ / _ \
    // | (_| | || (_) | | | | | | | | | | (_| | | | | | |  __/
    //  \__,_|\__\___/|_| |_| |_| |_| |_|\__,_|_| |_| |_|\___|
    */
    inline std::string Atom::atom_name() const
    {
        return _atom_name;
    }
    inline void Atom::set_atom_name(const std::string& s)
    {
        _atom_name = s;
    }

    /*        _ _     _
    //   __ _| | |_  | | ___   ___
    //  / _` | | __| | |/ _ \ / __|
    // | (_| | | |_  | | (_) | (__
    //  \__,_|_|\__| |_|\___/ \___|
    */
    inline char Atom::alt_loc() const
    {
        return _alt_loc;
    }
    inline void Atom::set_alt_loc(char a)
    {
        _alt_loc = a;
    }

    /*                _     _
    //  _ __ ___  ___(_) __| |  _ __   __ _ _ __ ___   ___
    // | '__/ _ \/ __| |/ _` | | '_ \ / _` | '_ ` _ \ / _ \
    // | | |  __/\__ \ | (_| | | | | | (_| | | | | | |  __/
    // |_|  \___||___/_|\__,_| |_| |_|\__,_|_| |_| |_|\___|
    */
    inline std::string Atom::resid_name() const
    {
        return _resid_name;
    }
    inline void Atom::set_resid_name(const std::string& s)
    {
        _resid_name = s;
    }

    /*      _           _         ___ ____
    //  ___| |__   __ _(_)_ __   |_ _|  _ \
    // / __| '_ \ / _` | | '_ \   | || | | |
    //| (__| | | | (_| | | | | |  | || |_| |
    // \___|_| |_|\__,_|_|_| |_| |___|____/
    */
    inline char Atom::chain_ID() const
    {
        return _chain_ID;
    }
    inline void Atom::set_chain_ID(char a)
    {
        _chain_ID = a;
    }

    /*                _     _   _           _
    //  _ __ ___  ___(_) __| | (_)_ __   __| | _____  __
    // | '__/ _ \/ __| |/ _` | | | '_ \ / _` |/ _ \ \/ /
    // | | |  __/\__ \ | (_| | | | | | | (_| |  __/>  <
    // |_|  \___||___/_|\__,_| |_|_| |_|\__,_|\___/_/\_\
    */
    inline unsigned int Atom::resid_index() const
    {
        return _resid_index;
    }
    inline void Atom::set_resid_index(unsigned int i)
    {
        _resid_index = i;
    }

    /*  _               _
    // (_) ___ ___   __| | ___
    // | |/ __/ _ \ / _` |/ _ \
    // | | (_| (_) | (_| |  __/
    // |_|\___\___/ \__,_|\___|
    */
    inline char Atom::icode() const
    {
        return _insert_code;
    }
    inline void Atom::set_icode(char a)
    {
        _insert_code = a;
    }

    /*                          _ _             _
    //   ___ ___   ___  _ __ __| (_)_ __   __ _| |_ ___  ___
    //  / __/ _ \ / _ \| '__/ _` | | '_ \ / _` | __/ _ \/ __|
    // | (_| (_) | (_) | | | (_| | | | | | (_| | ||  __/\__ \
    //  \___\___/ \___/|_|  \__,_|_|_| |_|\__,_|\__\___||___/
    */
    inline const Vec3d& Atom::coordinates() const
    {
        return _coordinate;
    }
    inline void Atom::set_coords(const Vec3d& coors)
    {
        _coordinate = coors;
    }

    /*            _            _ _   _
    // __   _____| | ___   ___(_) |_(_) ___  ___
    // \ \ / / _ \ |/ _ \ / __| | __| |/ _ \/ __|
    //  \ V /  __/ | (_) | (__| | |_| |  __/\__ \
    //   \_/ \___|_|\___/ \___|_|\__|_|\___||___/
    */
    inline const Vec3d& Atom::velocities() const
    {
        return _velocity;
    }
    inline void Atom::set_velocities(const Vec3d& velos)
    {
        _velocity = velos;
    }

    /*   ___   ___ ___ _   _ _ __   __ _ _ __   ___ _   _
    //  / _ \ / __/ __| | | | '_ \ / _` | '_ \ / __| | | |
    // | (_) | (_| (__| |_| | |_) | (_| | | | | (__| |_| |
    //  \___/ \___\___|\__,_| .__/ \__,_|_| |_|\___|\__, |
    //                      |_|                     |___/
    */
    inline double Atom::occupancy() const
    {
        return _occupancy;
    }
    inline void Atom::set_occupancy(double o)
    {
        _occupancy = o;
    }

    /*  _                          __            _
    // | |_ ___ _ __ ___  _ __    / _| __ _  ___| |_ ___  _ __
    // | __/ _ \ '_ ` _ \| '_ \  | |_ / _` |/ __| __/ _ \| '__|
    // | ||  __/ | | | | | |_) | |  _| (_| | (__| || (_) | |
    //  \__\___|_| |_| |_| .__/  |_|  \__,_|\___|\__\___/|_|
    //                   |_|
    */
    inline double Atom::temperature_factor() const
    {
        return _temp_factor;
    }
    inline void Atom::set_temperature_factor(double f)
    {
        _temp_factor = f;
    }

    /*                                      _     ___ ____
    //  ___  ___  __ _ _ __ ___   ___ _ __ | |_  |_ _|  _ \
    // / __|/ _ \/ _` | '_ ` _ \ / _ \ '_ \| __|  | || | | |
    // \__ \  __/ (_| | | | | | |  __/ | | | |_   | || |_| |
    // |___/\___|\__, |_| |_| |_|\___|_| |_|\__| |___|____/
    //           |___/
    */
    inline std::string Atom::segment_ID() const
    {
        return _seg_ID;
    }
    inline void Atom::set_segment_ID(const std::string& s)
    {
        _seg_ID = s;
    }

    /*       _                           _
    //   ___| | ___ _ __ ___   ___ _ __ | |_
    //  / _ \ |/ _ \ '_ ` _ \ / _ \ '_ \| __|
    // |  __/ |  __/ | | | | |  __/ | | | |_
    //  \___|_|\___|_| |_| |_|\___|_| |_|\__|
    */
    inline std::string Atom::element() const
    {
        return _element;
    }
    inline void Atom::set_element(const std::string& s)
    {
        _element = s;
    }

    /*       _
    //   ___| |__   __ _ _ __ __ _  ___
    //  / __| '_ \ / _` | '__/ _` |/ _ \
    // | (__| | | | (_| | | | (_| |  __/
    //  \___|_| |_|\__,_|_|  \__, |\___|
    //                       |___/
    */
    inline std::string Atom::charge() const
    {
        return _charge;
    }
    inline void Atom::set_charge(const std::string& s)
    {
        _charge = s;
    }

    // Atom::Atom ==============================================================
    inline Atom::Atom()
    {
        _atom_flag = "";
        _serial = 0;
        _atom_name = "";
        _alt_loc = ' ';
        _resid_name = "";
        _chain_ID = ' ';
        _resid_index = 0;
        _insert_code = ' ';
        _coordinate = Vec3d(0, 0, 0);
        _velocity = Vec3d(0, 0, 0);
        _occupancy = 0;
        _temp_factor = 0;
        _seg_ID = "";
        _element = "";
        _charge = "";
    }

    inline void Atom::reset()
    {
        _atom_flag = "";
        _serial = 0;
        _atom_name = "";
        _alt_loc = ' ';
        _resid_name = "";
        _chain_ID = ' ';
        _resid_index = 0;
        _insert_code = ' ';
        _coordinate = Vec3d(0, 0, 0);
        _velocity = Vec3d(0, 0, 0);
        _occupancy = 0;
        _temp_factor = 0;
        _seg_ID = "";
        _element = "";
        _charge = "";
    }

    /*   ___        _              _____                 _   _
    //  / _ \ _   _| |_ ___ _ __  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
    // | | | | | | | __/ _ \ '__| | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
    // | |_| | |_| | ||  __/ |    |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
    //  \___/ \__,_|\__\___|_|    |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
    */
    inline std::ostream& operator<<(std::ostream& o, const Atom& a)
    {
        if (a.atom_flag() == "ATOM  " || a.atom_flag() == "HETATM")
        {
            o << std::setw(6) << a.atom_flag()
              << std::setw(5) << a.serial() << " "
              << std::setw(4) << a.atom_name()
              << std::setw(1) << a.alt_loc()
              << std::setw(3) << a.resid_name() << " "
              << std::setw(1) << a.chain_ID()
              << std::setw(4) << a.resid_index()
              << std::setw(1) << a.icode() << "   "
              << a.coordinates()
              << std::setiosflags(std::ios_base::fixed) << std::setprecision(2)
              << std::setw(6) << a.occupancy()
              << std::setw(6) << a.temperature_factor() << "      "
              << std::setw(4) << a.segment_ID()
              << std::setw(2) << a.element()
              << std::setw(2) << a.charge();
        } else {
            o << a.atom_flag();
        }
        return o;
    }

    std::istream& operator>>(std::istream& i, Atom& a)
    {
        std::istringstream tmp_sstr;
        std::string pdb_line;

        unsigned int _tmp_ui;
        std::string _tmp_str;
        char _tmp_char;
        pinang::Vec3d _tmp_coordinates;
        double _tmp_d;

        std::getline( i, pdb_line);
        pdb_line.resize(80, ' ');
        _tmp_str = pdb_line.substr(0,6);
        a.set_atom_flag(_tmp_str);

        // a.set_atom_name("NULL");
        if (a.atom_flag() == "ATOM  " || a.atom_flag() == "HETATM")
        {
            tmp_sstr.str ( pdb_line.substr(6,5));
            tmp_sstr >> _tmp_ui;
            a.set_serial(_tmp_ui);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(12,4));
            tmp_sstr >> _tmp_str;
            a.set_atom_name(_tmp_str);
            tmp_sstr.clear();
            _tmp_str.clear();

            tmp_sstr.str ( pdb_line.substr(16,1));
            tmp_sstr.get(_tmp_char);
            a.set_alt_loc(_tmp_char);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(17,3));
            tmp_sstr >> _tmp_str;
            a.set_resid_name(_tmp_str);
            tmp_sstr.clear();
            _tmp_str.clear();

            tmp_sstr.str ( pdb_line.substr(21,1));
            tmp_sstr >> _tmp_char;
            a.set_chain_ID(_tmp_char);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(22,4));
            tmp_sstr >> _tmp_ui;
            a.set_resid_index(_tmp_ui);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(26,1));
            tmp_sstr.get(_tmp_char);
            a.set_icode(_tmp_char);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(30,24));
            tmp_sstr >> _tmp_coordinates;
            a.set_coords(_tmp_coordinates);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(54,6));
            tmp_sstr >> _tmp_d;
            a.set_occupancy(_tmp_d);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(60,6));
            tmp_sstr >> _tmp_d;
            a.set_temperature_factor(_tmp_d);
            tmp_sstr.clear();

            tmp_sstr.str ( pdb_line.substr(72,4));
            tmp_sstr >> _tmp_str;
            a.set_segment_ID(_tmp_str);
            tmp_sstr.clear();
            _tmp_str.clear();

            tmp_sstr.str ( pdb_line.substr(76,2));
            tmp_sstr >> _tmp_str;
            a.set_element(_tmp_str);
            tmp_sstr.clear();
            _tmp_str.clear();

            tmp_sstr.str ( pdb_line.substr(78,2));
            tmp_sstr >> _tmp_str;
            a.set_charge(_tmp_str);
            tmp_sstr.clear();
            _tmp_str.clear();
        } else {
            // std::cerr << "ERROR: Wrong format!" << std::endl;
        }

        if (!i) return i;
        return i;
    }

    inline double atom_distance (Atom& a1, Atom& a2)
    {
        Vec3d v3 = a1.coordinates() - a2.coordinates();
        return v3.norm();
    }

}
#endif

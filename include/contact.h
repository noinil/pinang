// -*-c++-*-

#ifndef PINANG_CONTACT_H_
#define PINANG_CONTACT_H_

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "constants.h"

namespace pinang{
    class Contact
    {
    public:
        Contact(): _res_id_1(0), _res_id_2(0), _distance(0), _flag(0) {}
        Contact(int a, int b, double c);
        virtual ~Contact() {}

        inline void set(int a, int b, double c);

        inline int is_contact() const;

    protected:
        int _res_id_1, _res_id_2;
        double _distance;
        int _flag;
    };

    inline Contact::Contact(int a, int b, double c)
    {
        _res_id_1 = a;
        _res_id_2 = b;
        _distance = c;
        if (_distance < cutoff)
        {
            _flag = 1;
        } else {
            _flag = 0;
        }
    }

    inline void Contact::set(int a, int b, double c)
    {
        _res_id_1 = a;
        _res_id_2 = b;
        _distance = c;
        if (_distance < cutoff)
        {
            _flag = 1;
        } else {
            _flag = 0;
        }
    }

    inline int Contact::is_contact() const
    {
        if (_res_id_1 * _res_id_2 == 0)
        {
            std::cerr << "ERROR: This contact is not properly initialized!"
                 << std::endl;
            exit(EXIT_SUCCESS);
        } else {
            return _flag;
        }
    }

}

#endif

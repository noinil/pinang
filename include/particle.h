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

  inline std::string get_atom_name() const;
  inline void set_atom_name(const std::string& s);

  inline std::string get_resid_name() const;
  inline void set_resid_name(const std::string& s);

  inline int get_resid_index() const;
  inline void set_resid_index(int i);

  inline double get_charge() const;
  inline void set_charge(double c);

  inline double get_mass() const;
  inline void set_mass(double m);

 protected:
  std::string atom_name_;
  std::string resid_name_;
  int resid_index_;
  double charge_;
  double mass_;
};

/*        _
//   __ _| |_ ___  _ __ ___    _ __   __ _ _ __ ___   ___
//  / _` | __/ _ \| '_ ` _ \  | '_ \ / _` | '_ ` _ \ / _ \
// | (_| | || (_) | | | | | | | | | | (_| | | | | | |  __/
//  \__,_|\__\___/|_| |_| |_| |_| |_|\__,_|_| |_| |_|\___|
*/
inline std::string Particle::get_atom_name() const
{
  return atom_name_;
}
inline void Particle::set_atom_name(const std::string& s)
{
  size_t sz = 3;
  atom_name_ = s;

  if (atom_name_.size() < sz)
  {
    atom_name_.resize(sz, ' ');
  }
}

/*                _     _
//  _ __ ___  ___(_) __| |  _ __   __ _ _ __ ___   ___
// | '__/ _ \/ __| |/ _` | | '_ \ / _` | '_ ` _ \ / _ \
// | | |  __/\__ \ | (_| | | | | | (_| | | | | | |  __/
// |_|  \___||___/_|\__,_| |_| |_|\__,_|_| |_| |_|\___|
*/
inline std::string Particle::get_resid_name() const
{
  return resid_name_;
}
inline void Particle::set_resid_name(const std::string& s)
{
  resid_name_ = s;
}

/*                _     _   _           _
//  _ __ ___  ___(_) __| | (_)_ __   __| | _____  __
// | '__/ _ \/ __| |/ _` | | | '_ \ / _` |/ _ \ \/ /
// | | |  __/\__ \ | (_| | | | | | | (_| |  __/>  <
// |_|  \___||___/_|\__,_| |_|_| |_|\__,_|\___/_/\_\
*/
inline int Particle::get_resid_index() const
{
  return resid_index_;
}
inline void Particle::set_resid_index(int i)
{
  resid_index_ = i;
}

/*       _
//   ___| |__   __ _ _ __ __ _  ___
//  / __| '_ \ / _` | '__/ _` |/ _ \
// | (__| | | | (_| | | | (_| |  __/
//  \___|_| |_|\__,_|_|  \__, |\___|
//                       |___/
*/
inline double Particle::get_charge() const
{
  return charge_;
}
inline void Particle::set_charge(double c)
{
  charge_ = c;
}

inline double Particle::get_mass() const
{
  return mass_;
}
inline void Particle::set_mass(double m)
{
  mass_ = m;
}

// Particle::Particle ==============================================================
inline Particle::Particle()
{
  atom_name_ = "";
  resid_name_ = "";
  resid_index_ = 0;
  charge_ = 0;
  mass_ = 0;
}

inline void Particle::reset()
{
  atom_name_ = "";
  resid_name_ = "";
  resid_index_ = 0;
  charge_ = 0;
  mass_ = 0;
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
  // std::cout << ".";

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

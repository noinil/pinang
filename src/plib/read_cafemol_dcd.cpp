#include "read_cafemol_dcd.hpp"

namespace pinang {

int read_cafemol_dcd(std::ifstream& dcd_file, std::vector<Conformation>& cfms)
{
  const std::size_t Si = sizeof(int);
  const std::size_t Sf = sizeof(float);
  // const std::size_t Sd = sizeof(double);
  // const std::size_t Sc = sizeof(char);

  int filesize = 0;
  int i_tmp = 0;
  float f_tmp = 0;

  std::vector<Vec3d> vv_tmp;
  Vec3d v_tmp(0,0,0);

  // ---------------------------------------------------------------------
  // DCD variables;
  // ---------------------------------------------------------------------
  // block 1;
  int flag1 = 0;
  int flag2 = 0;
  int flag3 = 0;
  int flag4 = 0;
  char * hdr_buf = new char [4];
  int nset = 0;
  int istart = 0;
  int nsavc = 0;
  int nstep = 0;
  int nunit = 0;
  int nfreat = 0;
  float delta = 0;
  int xtc = 0;
  int nver = 0;

  // block 2;
  int ntitle = 0;

  // block 3;
  int natom = 0;

  /*   __ _ _             _               _
  //  / _(_) | ___    ___| |__   ___  ___| | __
  // | |_| | |/ _ \  / __| '_ \ / _ \/ __| |/ /
  // |  _| | |  __/ | (__| | | |  __/ (__|   <
  // |_| |_|_|\___|  \___|_| |_|\___|\___|_|\_\
  */
  if (!dcd_file.is_open())
  {
    std::cout << " ~               PINANG :: DCD                ~ " << std::endl;
    std::cerr << " ERROR: Cannot read dcd file: not open." << std::endl;
    std::cerr << " Program terminating." << std::endl;
    return 1;
  }

  dcd_file.seekg(0, dcd_file.end);
  filesize = dcd_file.tellg();
  if (filesize == 0)
  {
    std::cout << " ~               PINANG :: DCD                ~ " << std::endl;
    std::cout << " !!! Error: empty dcd file! !!!" << std::endl;
    return 1;
  }
  dcd_file.seekg(0, dcd_file.beg);

  // ---------------------------------------------------------------------
  /*                     _   _     _            _      _
  //  _ __ ___  __ _  __| | | |__ | | ___   ___| | __ / |
  // | '__/ _ \/ _` |/ _` | | '_ \| |/ _ \ / __| |/ / | |
  // | | |  __/ (_| | (_| | | |_) | | (_) | (__|   <  | |
  // |_|  \___|\__,_|\__,_| |_.__/|_|\___/ \___|_|\_\ |_|
  */
  // ---------------------------------------------------------------------
  // First thing in the file should be an 84
  dcd_file.read((char *)&flag1, Si);
  if (flag1 != 84)
  {
    std::cout << " ~               PINANG :: DCD                ~ " << std::endl;
    std::cout << " !!! Magic number of block 1 error. !!!" << std::endl;
    return 1;
  }
  // ---------------------------------------------------------------------
  // 'CORD' for coordinate, 'VELD' for velocity
  dcd_file.read(hdr_buf ,4);
  // ---------------------------------------------------------------------
  // Store the number of sets of coordinates (NSET)
  dcd_file.read((char*)&nset, Si);
  // ---------------------------------------------------------------------
  // Store ISTART, the starting timestep
  dcd_file.read((char*)&istart, Si);
  // ---------------------------------------------------------------------
  // Store NSAVC, the number of timesteps between dcd saves
  dcd_file.read((char*)&nsavc, Si);
  // ---------------------------------------------------------------------
  // Store NSTEP, the number of steps
  dcd_file.read((char*)&nstep, Si);
  // ---------------------------------------------------------------------
  // Store NUNIT, the number of unit
  dcd_file.read((char*)&nunit, Si);
  // ---------------------------------------------------------------------
  // null
  for (int i = 0; i < 3; i++) {
    dcd_file.read((char*)&i_tmp, Si);
  }
  // ---------------------------------------------------------------------
  // Store NFREAT, the number of free atoms
  dcd_file.read((char*)&nfreat, Si);
  // ---------------------------------------------------------------------
  // Store DELTA, time step
  dcd_file.read((char*)&delta, Sf);
  // ---------------------------------------------------------------------
  // Store unit-cell info
  dcd_file.read((char*)&xtc, Si);
  // ---------------------------------------------------------------------
  // null
  for (int i = 0; i < 8; i++) {
    dcd_file.read((char*)&i_tmp, Si);
  }
  // Store NVER, version of CHARMM
  dcd_file.read((char*)&nver, Si);
  // ---------------------------------------------------------------------
  // Block size of the first block
  dcd_file.read((char*)&flag1, Si);
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  /*                     _   _     _            _      ____
  //  _ __ ___  __ _  __| | | |__ | | ___   ___| | __ |___ \
  // | '__/ _ \/ _` |/ _` | | '_ \| |/ _ \ / __| |/ /   __) |
  // | | |  __/ (_| | (_| | | |_) | | (_) | (__|   <   / __/
  // |_|  \___|\__,_|\__,_| |_.__/|_|\___/ \___|_|\_\ |_____|
  */
  // ---------------------------------------------------------------------
  // Block size of the second block
  dcd_file.read((char*)&flag2, Si);
  // ---------------------------------------------------------------------
  // Store NTITLE, the line number of title lines
  dcd_file.read((char*)&ntitle, Si);
  // ---------------------------------------------------------------------
  // Read in title
  for (int i=0; i<ntitle; i++)
    dcd_file.seekg(80, dcd_file.cur);
  // ---------------------------------------------------------------------
  // Block size of the second block
  dcd_file.read((char*)&flag2, Si);
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  /*                     _   _     _            _      _____
  //  _ __ ___  __ _  __| | | |__ | | ___   ___| | __ |___ /
  // | '__/ _ \/ _` |/ _` | | '_ \| |/ _ \ / __| |/ /   |_ \
  // | | |  __/ (_| | (_| | | |_) | | (_) | (__|   <   ___) |
  // |_|  \___|\__,_|\__,_| |_.__/|_|\___/ \___|_|\_\ |____/
  */
  // ---------------------------------------------------------------------
  // Block size of the third block; should be  4 !!!
  dcd_file.read((char*)&flag3, Si);
  if (flag3 != 4)
  {
    std::cout << " ~               PINANG :: DCD                ~ " << std::endl;
    std::cout << " !!! Magic number of block 3 error. !!!" << std::endl;
    return 1;
  }
  // ---------------------------------------------------------------------
  // Store NATOM, the number of atoms
  dcd_file.read((char*)&natom, Si);
  vv_tmp.clear();
  for (int i = 0; i < natom; ++i) {
    vv_tmp.push_back(v_tmp); // set the temp vector_of_vec3d; ~~~~~~
  }
  Conformation cfm_tmp(vv_tmp);
  // ---------------------------------------------------------------------
  // Block size of the third block; should be  4 !!!
  dcd_file.read((char*)&flag3, Si);
  if (flag3 != 4)
  {
    std::cout << " ~               PINANG :: DCD                ~ " << std::endl;
    std::cout << " !!! Magic number of block 3 error. !!!" << std::endl;
    return 1;
  }
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  /*                     _    ____ ___   ___  ____
  //  _ __ ___  __ _  __| |  / ___/ _ \ / _ \|  _ \
  // | '__/ _ \/ _` |/ _` | | |  | | | | | | | |_) |
  // | | |  __/ (_| | (_| | | |__| |_| | |_| |  _ <
  // |_|  \___|\__,_|\__,_|  \____\___/ \___/|_| \_\
  */
  // ---------------------------------------------------------------------
  while (!dcd_file.eof()) {
    dcd_file.read((char*)&flag4, Si);
    if (flag4/4 != natom)
    {
      std::cout << " !!! Coordinates mismatch the atom number. !!! "
                << std::endl;
      return 1;
    }
    // Read in all of X:
    for (int i = 0; i < natom; ++i) {
      dcd_file.read((char*)&f_tmp, Sf);
      vv_tmp[i] += Vec3d(f_tmp, 0, 0);
    }
    dcd_file.read((char*)&flag4, Si);
    // Read in all of Y:
    dcd_file.read((char*)&flag4, Si);
    for (int i = 0; i < natom; ++i) {
      dcd_file.read((char*)&f_tmp, Sf);
      vv_tmp[i] += Vec3d(0, f_tmp, 0);
    }
    dcd_file.read((char*)&flag4, Si);
    // Read in all of Z:
    dcd_file.read((char*)&flag4, Si);
    for (int i = 0; i < natom; ++i) {
      dcd_file.read((char*)&f_tmp, Sf);
      vv_tmp[i] += Vec3d(0, 0, f_tmp);
    }
    dcd_file.read((char*)&flag4, Si);

    if (dcd_file.eof())
      break;

    cfm_tmp.set_conformation(vv_tmp); // --- set temp conformation ---
    cfms.push_back(cfm_tmp);

    // ---------- reset temp vectors ----------
    vv_tmp.clear();
    v_tmp = v_tmp * 0;
    for (int i = 0; i < natom; ++i) {
      vv_tmp.push_back(v_tmp); // set the temp vector_of_vec3d; ~~~~~~
    }

  }
  // std::cout << cfms.size() << std::endl;

  delete [] hdr_buf;

  return 0;
}

}  // pinang

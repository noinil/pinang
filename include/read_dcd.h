// -*-c++-*-

#ifndef PINANG_READ_DCD_H_
#define PINANG_READ_DCD_H_

#include "constants.h"
#include "vec3d.h"

#include <vector>
#include <fstream>

namespace pinang {
    int read_dcd(std::ifstream& dcd_file, std::vector<Vec3d>& v3d)
    {
        const std::size_t Si = sizeof(int);
        const std::size_t Sd = sizeof(double);
        int filesize;
        char * str_tmp;

        dcd_file.seekg(0, dcd_file.end);
        filesize = dcd_file.tellg();
        if (filesize == 0)
        {
            std::cout << " *** Error: empty dcd file! ***" << std::endl;
            return 1;
        }
        dcd_file.seekg(0, dcd_file.beg);

        str_tmp = new char[Si];
        dcd_file.read(str_tmp, sizeof(int));
        int flag1 = std::atoi(str_tmp);
        std::cout << flag1 << std::endl;
    }

}
#endif

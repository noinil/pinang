#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
              << std::endl;
    std::cout << " ~           PINANG cafemol .ts file re-output            ~ "
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;
    std::cout << " Usage: "
              << argv[0]
              << " -f some.ts [-o some.dat] [-h]"
              << std::endl;
    std::cout << " ========================================================== "
              << std::endl;

    int opt;
    int ts_flag = 0;

    std::string ts_name = "some.ts";
    std::string out_name = "some.dat";

    while ((opt = getopt(argc, argv, "f:o:h")) != -1) {
        switch (opt) {
        case 'f':
            ts_name = optarg;
            ts_flag = 1;
            break;
        case 'o':
            out_name = optarg;
            break;
        case 'h':
            std::cout << " This program - "
                      << argv[0] << " - "
                      << "is used to trim the output of .ts generated by CafeMol."
                      << std::endl;
            exit(EXIT_SUCCESS);
            break;
        default: /* '?' */
            std::cout << " This program - "
                      << argv[0] << " - "
                      << "is used to trim the output of .ts generated by CafeMol."
                      << std::endl
                      << " Please provide the .ts file. "
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (ts_flag == 0)
    {
        std::cout << " ERROR: Please provide the .ts file. " << std::endl;
        exit(EXIT_SUCCESS);
    }
    // -------------------------------------------------------------------------

    std::ifstream ts_file(ts_name.c_str());
    std::ofstream out_file(out_name.c_str());

    /*                  _         _
    //  _ __ ___   __ _(_)_ __   | | ___   ___  _ __
    // | '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \
    // | | | | | | (_| | | | | | | | (_) | (_) | |_) |
    // |_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/
    //                                         |_|
    */

    std::string inp_line;
    std::istringstream tmp_str_str;

    std::string str0 = "#unit";
    std::string str1 = "#all";
    std::string str_tmp = "#";
    char c_tmp = ' ';

    std::cout << " Please choose your molecule: \n"
              << "  - 0 : #all\n"
              << "  - 1 : #1\n"
              << "  - 2 : #2\n"
              << "  - ... "
              << std::endl << " ";
    std::cin >> c_tmp;
    if (c_tmp > '0')
    {
        if (c_tmp > '9')
        {
            std::cout << " ERROR: Please input an inteer. " << std::endl;
            return 1;
        }
        str1 = "#";
        str1.push_back(c_tmp);
        std::cout << " You choosed molecule " << str1
                  << ", Master!"
                  << std::endl
                  << "  ......"
                  << std::endl;
    }

    int i_count = 0;
    int i_mol = 0;
    std::size_t i_found;
    std::vector < int > v_quantity_num;
    while (ts_file.good()) {
        std::getline(ts_file, inp_line);
        if (ts_file.fail())
            break;

        // ==================== read in the #unit line ====================
        i_found = inp_line.find(str0);
        if (i_found!=std::string::npos){
            std::cout << " Units found!  "
                      << "Please choose your quantities: (end with 0)"
                      << std::endl;
            tmp_str_str.str(inp_line);
            std::vector<std::string> v_str;
            while (tmp_str_str.good()) {
                tmp_str_str >> str_tmp;
                if (tmp_str_str.fail())
                    break;
                if (str_tmp == "#unit" || str_tmp == "step")
                    continue;
                i_count++;
                v_str.push_back(str_tmp);
                std::cout << " " << i_count << " : "
                          << str_tmp << std::endl;
            }
            tmp_str_str.clear();

            std::cout << " " ;
            int ii = 0;
            while (cin >> ii) {
                if (ii == 0)
                    break;
                if ( ii > i_count) {
                    std::cout << " ERROR: Quantity No. out of range! "
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                v_quantity_num.push_back(ii);
            }
            std::cout << " Now I'm trying to extract these quantities: ";
            for (unsigned int i = 0; i < v_quantity_num.size(); i++)
                std::cout << v_str[v_quantity_num[i]-1] << " ";
            std::cout << std::endl;

            out_file << "#" << std::setw(9) << "step";
            for (unsigned int i = 0; i < v_quantity_num.size(); i++)
                out_file << std::setw(12) << v_str[v_quantity_num[i]-1];
            out_file << std::endl;
        }

        i_found = inp_line.find(str1);
        if (i_found!=std::string::npos){
            i_mol = 1;
            long int time_step;
            if (i_count == 0)
            {
                std::cout << " ERROR: No #unit line found before quantities."
                          << std::endl;
                exit(EXIT_FAILURE);
            }
            tmp_str_str.str(inp_line);
            tmp_str_str >> str_tmp; // read in "#all" or "#1" ...;
            tmp_str_str >> time_step; // read in step;

            double dd = 0;      // temp read-in double;
            std::vector<double> vd; // tmp vector for quantities;
            while (tmp_str_str.good()) {
                if (tmp_str_str.fail())
                    break;
                tmp_str_str >> dd;
                vd.push_back(dd);
            }

            out_file << std::setw(10) << time_step;
            for (unsigned int i = 0; i < v_quantity_num.size(); i++)
                out_file << std::setw(12) << vd[v_quantity_num[i]-1];
            out_file << std::endl;
            vd.clear();
            tmp_str_str.clear();
        }
    }
    if (i_mol == 0)
    {
        std::cout << " ERROR: You have chosen the Wrong Molecule! " << std::endl;
    }

    ts_file.close();
    out_file.close();

    return 0;
}

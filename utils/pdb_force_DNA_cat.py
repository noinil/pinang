#!/usr/bin/env python

import sys

def main(p0_name):
    """Change single char names of DNA atoms to "D*" style.
    p0_name: original PDB file name.
    output:  a new PDB
    """
    p1_name = p0_name[:-4] + '_DNAnization.pdb'
    p1 = open(p1_name, 'w')

    keyword_list = ['  A', '  G', '  C', '  T', \
                    'A  ', 'G  ', 'C  ', 'T  ', \
                    ' A ', ' G ', ' C ', ' T ']
    new_name_list = [' DA', ' DG', ' DC', ' DT']

    with open(p0_name, 'r') as fin:
        for line in fin:
            keyword = line[17:20]
            if keyword in keyword_list:
                new_name = new_name_list[keyword_list.index(keyword) % 4]
            else:
                new_name = keyword
            new_line = line[:17] + new_name + line[20:]
            p1.write(new_line)

    p1.close()


if __name__ == '__main__':
    pdb_0_name = sys.argv[1]
    main(pdb_0_name)

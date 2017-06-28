#!/usr/bin/env python3

def main(cafe_pdb_fname):
    fout_name = cafe_pdb_fname + '.normal.pdb'
    fout = open(fout_name, 'w')
    with open(cafe_pdb_fname, 'r') as fin:
        final_shot_flag = False
        for line in fin:
            if final_shot_flag:
                if line.startswith('>>'):
                    fout.write('TER\n')
                elif line.startswith('ATOM'):
                    fout.write(line)
                elif line.startswith('ENDMDL'):
                    fout.write('END\n')
            if line.startswith("MODEL"):
                if line.split()[1] == '2':
                    final_shot_flag = True
    fout.close()
    

if __name__ == '__main__':
    import sys
    main(sys.argv[1])

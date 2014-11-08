#!/usr/bin/env python

n_bp = 0

def main():
    global n_bp

    bond_dict = {}
    angl_dict = {}
    dihe_dict = {}

    f_in_bond = open("in00_bond.list", 'r')
    f_in_angl = open("in00_angl.list", 'r')
    f_in_dihe = open("in00_dihe.list", 'r')
    f_strnd_1 = open("strand1.ninfo", 'w')
    f_strnd_2 = open("strand2.ninfo", 'w')

    f_strnd_1.write("<<<< native bond length\n")
    f_strnd_2.write("<<<< native bond length\n")

    bond_string = "bond   {:4d}     {:2d}     {:2d}   {:4d}   {:4d}    {:4d}    {:4d} \
    {:8.4f}     {:8.4f}     {:8.4f}     {:8.4f}\n"
    angl_string = "angl   {:4d}     {:2d}     {:2d}   {:4d}   {:4d}    {:4d}    {:4d}    {:4d}    {:4d} \
    {:8.4f}     {:8.4f}     {:8.4f}     {:8.4f}\n"
    dihe_string = "dihd   {:4d}     {:2d}     {:2d}   {:4d}   {:4d}    {:4d}    {:4d}    {:4d}    {:4d}\
    {:4d}    {:4d}     {:8.4f}     {:8.4f}     {:8.4f}     {:8.4f}     {:8.4f}\n"

    # -------------------- BOND --------------------
    for lines in f_in_bond:
        words = lines.split()
        imp1, imp2 = int(words[0]), int(words[1])
        bd_nat = float(words[2])
        coef_bd = float(words[3])
        bond_dict[(imp1, imp2)] = [bd_nat, coef_bd]

    i_count = 0
    bond_keys = sorted(bond_dict.keys())
    for k in bond_keys:
        i_count += 1
        v = bond_dict[k]
        if k[0] > n_bp or k[1] > n_bp:
            f_strnd_2.write(bond_string.format(i_count, 2, 2, k[0], k[1], k[0] - n_bp, \
                                               k[1] - n_bp, v[0], 1.0, 1.0, v[1]))
        else:
            f_strnd_1.write(bond_string.format(i_count, 1, 1, k[0], k[1], k[0], \
                                               k[1], v[0], 1.0, 1.0, v[1]))
    f_strnd_1.write(">>>>\n\n")
    f_strnd_2.write(">>>>\n\n")

    # -------------------- ANGLE --------------------
    f_strnd_1.write("<<<< native bond angles\n")
    f_strnd_2.write("<<<< native bond angles\n")
    for lines in f_in_angl:
        words = lines.split()
        imp1, imp2, imp3 = int(words[0]), int(words[1]), int(words[2])
        ba_nat = float(words[3])
        coef_ba = float(words[4])
        angl_dict[(imp1, imp2, imp3)] = [ba_nat, coef_ba]
    i_count = 0
    angl_keys = sorted(angl_dict.keys())
    for k in angl_keys:
        i_count += 1
        v = angl_dict[k]
        if k[0] > n_bp or k[1] > n_bp or k[2] > n_bp:
            f_strnd_2.write(angl_string.format(i_count, 2, 2, k[0], k[1], k[2], k[0] - n_bp, \
                                               k[1] - n_bp, k[2] - n_bp, v[0], 1.0, 1.0, v[1]))
        else:
            f_strnd_1.write(angl_string.format(i_count, 1, 1, k[0], k[1], k[2], k[0], \
                                               k[1], k[2], v[0], 1.0, 1.0, v[1]))
    f_strnd_1.write(">>>>\n\n")
    f_strnd_2.write(">>>>\n\n")

    # -------------------- ANGLE --------------------
    f_strnd_1.write("<<<< native dihedral angles\n")
    f_strnd_2.write("<<<< native dihedral angles\n")
    for lines in f_in_dihe:
        words = lines.split()
        imp1, imp2, imp3, imp4 = int(words[0]), int(words[1]), int(words[2]), int(words[3])
        dih_nat = float(words[4])
        coef_prd = float(words[5])
        coef_gau = float(words[6])
        coef_sig = float(words[7])
        coef_type = float(words[0])
        dihe_dict[(imp1, imp2, imp3, imp4)] = [dih_nat, coef_prd, coef_gau, coef_sig, coef_type]
    i_count = 0
    dihe_keys = sorted(dihe_dict.keys())
    for k in dihe_keys:
        i_count += 1
        v = dihe_dict[k]
        if k[0] > n_bp or k[1] > n_bp or k[2] > n_bp:
            f_strnd_2.write(dihe_string.format(i_count, 2, 2, k[0], k[1], k[2], k[3], k[0] - n_bp, \
                                               k[1] - n_bp, k[2] - n_bp, k[3] - n_bp, v[0], 1.0, 1.0, v[2], v[3]))
        else:
            f_strnd_1.write(dihe_string.format(i_count, 1, 1, k[0], k[1], k[2], k[3], k[0], \
                                               k[1], k[2], k[3], v[0], 1.0, 1.0, v[2], v[3]))
    f_strnd_1.write(">>>>\n\n")
    f_strnd_2.write(">>>>\n\n")


    f_in_bond.close()
    f_in_angl.close()
    f_in_dihe.close()

if __name__ == '__main__':
    import sys
    filename = sys.argv[1]
    with open(filename, 'r') as f:
        fline = f.readline()
        n_bp = int(fline) * 3 - 1

    main()

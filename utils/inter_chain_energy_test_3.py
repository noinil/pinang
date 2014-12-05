#!/usr/bin/env python

T = 300.0
I = 0.1
ek = 78.0                       # diele const
coef_exv = 0.1
coef_exv_hard = 0.2
B = 12
C = 6

def distance(xi, xj):
    x1, y1, z1 = xi
    x2, y2, z2 = xj
    return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5

def normalize(vec):
    x, y, z = vec[0], vec[1], vec[2]
    r = (x*x + y*y + z*z) ** 0.5
    return (x/r, y/r, z/r)

def pair_ene_ele(q1, q2, dist):
    import math
    global T
    global I
    global ek
    d_D = math.sqrt(3.95e-4*ek*T/I)
    # print(q1, q2, dist, d_D)
    E_ele = 322 * q1 * q2 * math.exp(-dist/ d_D) / (ek * dist)
    return E_ele

def pair_ene_exv(dist, sigma):
    import math
    global coef_exv
    global coef_exv_hard
    global B
    global C
    c_tmp = sigma/dist
    c_tmp2 = 4.0/dist
    if c_tmp <= 1:
        E_exv = 0
    else:
        E_exv = coef_exv * (c_tmp ** B - (B/C) * c_tmp ** C + (B/C - 1))
    # if c_tmp2 <= 1:
    #     E_exv += 0
    # else:
    # E_exv += coef_exv_hard * (c_tmp2**12)
    return E_exv

def main(top_name, pos_name):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab

    coors = []                  # coordinates of every resid
    charges = []                # charge of every residue
    table_ex_rad = []                 # excluded volumn radii
    mol_type = []               # molecular type
    chain_ind = []              # resid in which chain?
    table_sigma_dist = []
    min_dist_matrix = []

    for i in range(50):
        sig_dist = {}
        q = 0.01 * i + 0.01
        q_str = str(round(q,2))
        if len(q_str) < 4:
            q_str += '0'
        sigma_file = open("/home/noinil/Workspace/excluded_volumn_test/final_radius/sigma_out_"+q_str+".dat", 'r')
        for lines in sigma_file:
            str0, radi = lines.split()[0], float(lines.split()[1])
            if str0 == 'A':
                radi = 5.4
            elif str0 == 'T':
                radi = 7.1
            elif str0 == 'G':
                radi = 4.9
            elif str0 == 'C':
                radi = 6.4
            elif str0 == 'P':
                radi = 4.5
            elif str0 == 'S':
                radi = 6.2
            else:
                radi = radi
            sig_dist[str0] = radi
        table_sigma_dist.append(sig_dist)

    with open(top_name, 'r') as top_f:
        read_flag = 0
        for lines in top_f:
            words = lines.split()
            if len(words) > 1 and words[1] == "particles":
                read_flag = 1
            if len(words) > 1 and words[1] == "bonds":
                break
            if read_flag == 0:
                continue
            if len(words) == 6:
                if words[-1] == "-0.6":
                    charge = -1.0
                else:
                    charge = float(words[-1])
                charges.append(charge)
                if words[3] == "CA":
                    mol_type.append(0)
                else:
                    mol_type.append(1)

    first_pro_chain_ind = -1
    res_name = []
    with open(pos_name, 'r') as pos_f:
        chain_num = 0
        for lines in pos_f:
            words = lines.split()
            if len(words) > 2 and words[1] == "Chain":
                chain_num += 1
            if len(words) >= 6 and words[0].isnumeric():
                resnum = int(words[0])
                resname = words[1]
                x, y, z = float(words[3]), float(words[4]), float(words[5])
                res_name.append(resname)
                if mol_type[resnum - 1] == 0 and first_pro_chain_ind == -1:
                    first_pro_chain_ind = chain_num # ========== MARK the 1st pro chain
                chain_ind.append(chain_num)
                coors.append((x, y, z))
    for i in range(50):
        ex_rad = []
        for j in res_name:
            ex_rad.append(table_sigma_dist[i][j])
        table_ex_rad.append(ex_rad[:])
        ex_rad.clear()
    # for i, k in enumerate(table_sigma_dist):
    #     print(i, k)
    # for i in table_ex_rad:
    #     print(i)
    print(" Chain num =", chain_num)
    if chain_num > 3:
        print("  might be wrong, \'cause chain number is more than 3...")
        print("  calculating using only the first chain:", first_pro_chain_ind)

    # ==================== neighboring list ====================
    charge_list = []
    exv_list = []
    nearest_group = [set(), set()]
    for i, u in enumerate(mol_type):
        if u == 0 and chain_ind[i] != first_pro_chain_ind:
            continue
        for j, v in enumerate(mol_type):
            if v == 0 and chain_ind[j] != first_pro_chain_ind:
                continue
            if u != v and i < j:
                c1, c2 = charges[i], charges[j]
                dist = distance(coors[i], coors[j])
                if abs(c1 * c2) > 0.1:
                    charge_list.append((i, j, c1, c2))
                sig0 = (table_ex_rad[-1][i] + table_ex_rad[-1][j]) / 2
                if dist < sig0 + 10:
                    exv_list.append((i, j))
                    nearest_group[u].add(i)
                    nearest_group[v].add(j)

    # ---------- get the COM of pro and DNA interface ----------
    x, y, z = 0, 0, 0
    for m in nearest_group[0]:
        x += coors[m][0]
        y += coors[m][1]
        z += coors[m][2]
    x /= len(nearest_group[0])
    y /= len(nearest_group[0])
    z /= len(nearest_group[0])
    com0 = (x, y, z)            # ---------- COM of pro ----------

    x, y, z = 0, 0, 0
    for m in nearest_group[1]:
        x += coors[m][0]
        y += coors[m][1]
        z += coors[m][2]
    x /= len(nearest_group[1])
    y /= len(nearest_group[1])
    z /= len(nearest_group[1])
    com1 = (x, y, z)            # ---------- COM of DNA ----------

    vec_01 = (com1[0]-com0[0], com1[1]-com0[1], com1[2]-com0[2])
    vec_01 = normalize(vec_01)  # From vec0 to vec 1

    X = [-2 + 0.1 * i for i in range(60)]
    # table_E_TOT = []
    # for p in range(10):
    for p in range(20):
        global coef_exv
        coef_exv = p * 0.05 + 0.05
        # coef_exv = p * 0.02 + 0.02
        min_dist = []
        for q in range(20):
            # quantile = q * 5 + 4
            quantile = q
            ex_rad = table_ex_rad[quantile]
            E_TOT = []
            for i in range(30):
                dr = -2 + i * 0.2
                dx, dy, dz = dr * vec_01[0], dr * vec_01[1], dr * vec_01[2]
                tmp_coors = []
                for j, ctmp in enumerate(coors):
                    x0, y0, z0 = ctmp
                    if mol_type[j] == 0 and chain_ind[j] == first_pro_chain_ind:
                        tmp_coors.append((x0-dx, y0-dy, z0-dz))
                    else:
                        tmp_coors.append((x0, y0, z0))

                # -------------------- calc energies --------------------
                total_energy = 0
                for k in charge_list:
                    imp, jmp = k[0], k[1]
                    chi, chj = k[2], k[3]
                    dist = distance(tmp_coors[imp], tmp_coors[jmp])
                    e_tmp =  pair_ene_ele(chi, chj, dist)
                    total_energy += e_tmp
                for k in exv_list:
                    imp, jmp = k[0], k[1]
                    dist = distance(tmp_coors[imp], tmp_coors[jmp])
                    sig = 0.5 * (ex_rad[imp] + ex_rad[jmp])
                    e_tmp = pair_ene_exv(dist, sig)
                    total_energy += e_tmp
                tmp_coors.clear()
                E_TOT.append(total_energy)
            i_min, e_min = 10, E_TOT[0] + 1
            for i,e in enumerate(E_TOT):
                if e < e_min:
                    i_min, e_min = i, e
                else:
                    break
            E_TOT.clear()
            min_dist.append(i_min * 0.2 - 2)
        min_dist_matrix.append(min_dist[:])
        min_dist.clear()

    matrix_data_file = open(top_name[:-4]+"_matrix.dat", 'w')
    for i in min_dist_matrix:
        for j in i:
            matrix_data_file.write(str(round(j, 3))+'  ')
        matrix_data_file.write('\n')

    x = [i for i in range(0,21,2)]
    xtx = [round(i*0.02 + 0.01, 2) for i in range(11)]
    ytx = [round(i*0.1 + 0.05, 2) for i in range(11)]
    plt.xlabel(r'quantile')
    plt.ylabel(r'$\epsilon$ (kcal/mol)')
    plt.xticks(x, xtx)
    plt.yticks(x, ytx)
    plt.title(top_name[:-4]+' ('+r'$N_{chain}=$'+str(chain_num)+')')
    cax = plt.imshow(min_dist_matrix, cmap=plt.cm.bwr, interpolation='none', origin='lower', \
                       vmin=-2.0, vmax=2.0)
    cbar = plt.colorbar(cax, ticks=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
    cbar.ax.set_yticklabels([r'-2$\AA$',  r'', r'-1$\AA$', r'',r'0$\AA$', r'', r'1$\AA$', r'', r'2$\AA$'])
    plt.savefig(top_name[:-4]+'_min_dist_matrix.png', dpi=150)
    # plt.show()

if __name__ == '__main__':
    import sys
    top_file_name = sys.argv[1]
    pos_file_name = sys.argv[2]
    main(top_file_name, pos_file_name)

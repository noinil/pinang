#!/usr/bin/env python

huge_repo = {}

def anaf(fin_name):
    """Read from distance.dat files.
    Keyword Arguments:
    fin_name -- data file name
    """
    with open(fin_name, 'r') as fin:
        for lines in fin:
            line = lines.split()
            if len(line) == 4 and line[0] == "CG_PAIR":
                resname1, resname2 = line[1], line[2]
                if len(resname1) < len(resname2):
                    print(" Exceptions of NA-AA pair! ")
                    resname1, resname2 = resname2, resname1
                if resname1 < resname2 and len(resname1) == len(resname2) == 3:
                    resname1, resname2 = resname2, resname1
                key0, dist0 = (resname1, resname2), float(line[3])
                if dist0 < 2.5 or dist0 > 20.0:
                    # print(fin_name, "  ", dist0, "  ", key0)
                    continue
                if key0 not in huge_repo:
                    huge_repo[key0] = []
                huge_repo[key0].append(dist0)

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    import glob
    for filename in glob.glob(r'./*.dat'):
        # print(filename)
        anaf(filename)

    pro_names = ['ALA', 'ARG', 'GLU', 'GLY', 'PHE', 'ASP', 'PRO', 'THR', 'SER',\
                'LYS', 'HIS', 'ILE', 'CYS', 'VAL', 'LEU', 'MET', 'TYR', 'TRP',\
                'GLN', 'ASN']
    na_names = ['A', 'T', 'G', 'C', 'P', 'S']

    # -------------------- INTER-PROTEIN --------------------
    pcommand = input(" Produce pro-pro distance distribution fig?  ")
    if pcommand == 'y' or pcommand == 'yes':
        for i in pro_names:
            fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(10,10))
            out_name = "inter_pro_" + i + ".png"

            for j, k in enumerate(pro_names):
                reskey = (i, k) if i > k else (k, i)
                print(reskey)
                m, n = j // 4, j % 4 # sub fig index
                x = np.array(huge_repo[reskey])
                num_bins = 20
                # the histogram of data
                p, bins, patches = axes[m, n].hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
                # y is Cumulative P
                y_sum = sum(p[:])
                y = [sum(p[:i])/y_sum for i in range(len(p)+1)]
                axes[m, n].plot(bins, y, 'r--')
                axes[m, n].grid(axis='y', linestyle='-', alpha=0.3)
                axes[m, n].set_xlabel(r'distance $\AA$')
                axes[m, n].set_ylabel('Probability')
                axes[m, n].set_xlim(0, 15)
                axes[m, n].set_ylim(0, 0.6)
                axes[m, n].set_xticks(np.arange(0, 15, 3))
                axes[m, n].set_title(i + ' - ' + k, fontsize=10)
            fig.subplots_adjust(hspace=0.8)
            fig.subplots_adjust(wspace=0.8)
            fig.savefig(out_name, dpi=150)

    # -------------------- PROTEIN -- DNA --------------------
    pcommand = input(" Produce pro-DNA distance distribution fig?  ")
    if pcommand == 'y' or pcommand == 'yes':
        for i in na_names:
            fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(10,10))
            out_name = "DNA_pro_" + i + ".png"

            for j, k in enumerate(pro_names):
                reskey = (k, i)
                print(reskey)
                m, n = j // 4, j % 4 # sub fig index
                x = np.array(huge_repo[reskey])
                num_bins = 20
                # the histogram of data
                p, bins, patches = axes[m, n].hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
                # y is Cumulative P
                y_sum = sum(p[:])
                y = [sum(p[:i])/y_sum for i in range(len(p)+1)]
                axes[m, n].plot(bins, y, 'r--')
                axes[m, n].grid(axis='y', linestyle='-', alpha=0.3)
                axes[m, n].set_xlabel(r'distance $\AA$')
                axes[m, n].set_ylabel('Probability')
                axes[m, n].set_xlim(0, 15)
                axes[m, n].set_ylim(0, 0.6)
                axes[m, n].set_xticks(np.arange(0, 15, 3))
                axes[m, n].set_title(i + ' - ' + k, fontsize=10)
            fig.subplots_adjust(hspace=0.8)
            fig.subplots_adjust(wspace=0.8)
            fig.savefig(out_name, dpi=150)

if __name__ == '__main__':
    main()

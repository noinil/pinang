#!/usr/bin/env python

huge_repo = {}                  # ['resA', 'resB'] : [d1, d2, d3, ...]
p_repo = {}                     # ['resA', 'resB'] : [(d1, p1), (d2, p2), ...]

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
    import matplotlib.axis as axis
    import glob
    for filename in glob.glob(r'./*.dat'):
        # print(filename)
        anaf(filename)

    pro_names = ['ARG', 'LYS', 'ASP', 'GLU', 'HIS',\
                 'SER', 'THR', 'GLN', 'ASN', 'ALA',\
                 'VAL', 'LEU', 'ILE', 'MET', 'PHE',\
                 'TYR', 'TRP', 'PRO', 'CYS', 'GLY']
    na_names = ['A', 'T', 'G', 'C', 'P', 'S']
    pna_names = pro_names + na_names

    # -------------------- INTER-PROTEIN --------------------
    pcommand = input(" Produce pro-pro distance distribution fig?  ")
    if pcommand == 'y' or pcommand == 'yes':
        for i in pro_names:
            fig, axes = plt.subplots(nrows=5, ncols=4, figsize=(10,10))
            out_name = "inter_pro_" + i + ".png"
            print(" Plotting ", out_name)

            for j, k in enumerate(pro_names):
                reskey = (i, k) if i > k else (k, i)
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
            print(" Plotting ", out_name)

            for j, k in enumerate(pro_names):
                reskey = (k, i)
                m, n = j // 4, j % 4 # sub fig index
                x = np.array(huge_repo[reskey])
                num_bins = 20
                # the histogram of data
                p, bins, patches = axes[m, n].hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
                # y is Cumulative P
                y_sum = sum(p[:])
                y = [sum(p[:i])/y_sum for i in range(len(p)+1)]
                # local_list = []
                # for kb, ky in zip(bins, y):
                #     local_list.append( (kb, ky) )
                # p_repo[reskey] = local_list

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

    pcommand = input(" Produce pro-DNA minimal distance map? (q to quit): ")
    if pcommand == 'q':
        return 0
    elif pcommand == 'y' or pcommand == 'yes':
        for k, l in huge_repo.items():
            sl = sorted(l)
            sl_L, sl_S = max(sl), min(sl)
            delta = (sl_L - sl_S) / 50
            y = [0 for i in range(51)]
            for dist in l:
                n_sl = int((dist - sl_S) / delta)
                y[n_sl] += 1
            y_sum = sum(y)
            for i in range(51):
                y[i] /= y_sum
            local_dict = [(sl_S + i * delta, y[i]) for i in range(50)]
            p_repo[k] = local_dict
        while True:
            p_cutoff_s = input(" Please give me a cutoff for P (q to quit): ")
            if p_cutoff_s == 'q':
                break
            p_cutoff = float(p_cutoff_s)

            imatrix = []
            for i, ni in enumerate(pna_names):
                row = []
                for j, nj in enumerate(pna_names):
                    reskey = (ni, nj) if (ni > nj and len(ni) == len(nj)) or len(ni) > len(nj) else (nj, ni)
                    if reskey not in p_repo:
                        # print(" Index error!", reskey, "not in p_repo!")
                        row.append(0)
                    else:
                        for dist_p in p_repo[reskey]:
                            if dist_p[1] > p_cutoff:
                                # print(dist_p[0])
                                break
                        row.append(dist_p[0])
                        # print(' end...')
                imatrix.append(row[:])
                row.clear()
            print(imatrix)
            # fig, ax = plt.subplots()
            x = [i for i in range(len(pna_names))]
            plt.xticks(x, pna_names, rotation='vertical')
            plt.yticks(x, pna_names)
            ax = plt.gca()
            ax.xaxis.set_stick_position('top')
            plt.imshow(imatrix, cmap=plt.cm.gray, interpolation='none', origin='upper')
            plt.show()


if __name__ == '__main__':
    main()

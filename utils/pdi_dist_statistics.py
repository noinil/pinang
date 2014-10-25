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
    import matplotlib
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
        pass
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
            local_dict = [(sl_S + i * delta, sum(y[:i])) for i in range(50)]
            p_repo[k] = local_dict[:]

        p_cutoff_list = [0.04 + 0.01 * i for i in range(20)]
        it = 0
        while True:
            # p_cutoff_s = input(" Please give me a cutoff for P (q to quit): ")
            # if p_cutoff_s == 'q':
            #     break
            # p_cutoff = float(p_cutoff_s)

            if it >= 20:
                break
            p_cutoff = p_cutoff_list[it]
            p_cutoff_s = str(round(p_cutoff, 2))
            it += 1
            print(it, "   ", p_cutoff)

            # -------------------- compute the distance matrix -----------------
            imatrix = []
            for i, ni in enumerate(pna_names):
                row = []
                for j, nj in enumerate(pna_names):
                    reskey = (ni, nj) if (ni > nj and len(ni) == len(nj)) or len(ni) > len(nj) else (nj, ni)
                    if reskey not in p_repo:
                        row.append(0)
                    else:
                        for dist_p in p_repo[reskey]:
                            if dist_p[1] > p_cutoff:
                                break
                        row.append(dist_p[0])
                imatrix.append(row[:])
                row.clear()
            # print(imatrix)

            x = [i for i in range(len(pna_names))]
            # -------------------- least square - compute sigmas ---------------
            # ---------- calc from prot first, then from prot-DNA --------------
            sigma, sigma_i_sum = [0 for i in range(26)], [0 for i in range(26)]
            total_sigma_sum = 0
            for i in range(20):
                sigma_sum = sum(imatrix[i][:20])
                sigma_i_sum[i] = sigma_sum
                total_sigma_sum += sigma_sum
            total_sigma_sum /= 40
            for i in range(20, 26):
                sigma_i_sum[i] = sum(imatrix[i][:20])
            for i in range(26):
                sigma[i] = 0.05 * (sigma_i_sum[i] - total_sigma_sum) * 2
            # ---------- calc directly from prot-DNA --------------
            sigma2, sigma_i_sum2 = [0 for i in range(26)], [0 for i in range(26)]
            total_sigma_sum2 = 0
            for i in range(20, 26):
                sigma_sum2 = sum(imatrix[i][:20])
                total_sigma_sum2 += sigma_sum2
            for i in range(26):
                sigma_i_sum2[i] = sum(imatrix[i][:])
            ts1 = total_sigma_sum + 0.05 * (total_sigma_sum2 - 6 * total_sigma_sum)
            for i in range(20):
                sigma2[i] = 2 * (sigma_i_sum2[i] - ts1) / 26
            ts2 = sum(sigma2[:20]) / 2
            for i in range(20, 26):
                sigma2[i] = 0.05 * (sigma_i_sum2[i] - ts2) * 2
            plt.plot(x, sigma, 'r-', linewidth=2, label='pro-pro')
            plt.plot(x, sigma2, 'g-', linewidth=2, label='pro-DNA')
            plt.xticks(x, pna_names, rotation='vertical')
            plt.ylim(3, 11)
            plt.ylabel(r'Distance ($\AA$)')
            plt.title(r'$\sigma$ (Quantile = '+p_cutoff_s+' )')
            plt.grid(axis='x', linestyle='--', alpha=0.3)
            plt.grid(axis='y', linestyle='-', alpha=0.3)
            plt.legend(prop={'size': 16}, loc='upper left')
            plt.savefig("sigma_quantile_"+p_cutoff_s+".png", dpi=150)
            # plt.show()
            plt.clf()

            # ==================== plot the matrix! ====================
            plt.xticks(x, pna_names, rotation='vertical')
            plt.yticks(x, pna_names)
            cax = plt.imshow(imatrix, cmap=plt.cm.BuGn, interpolation='none', origin='lower', \
                       vmin=2.0, vmax=10.0)
            plt.title("Distance Matrix Quantile = "+p_cutoff_s)
            cbar = plt.colorbar(cax, ticks=[2,4,6,8,10])
            cbar.ax.set_yticklabels([r'<2$\AA$', r'4$\AA$', r'6$\AA$', r'8$\AA$', r'10$\AA$'])
            plt.savefig("dist_matrix_quantile_"+p_cutoff_s+".png", dpi=150)
            # plt.show()
            plt.clf()


if __name__ == '__main__':
    main()

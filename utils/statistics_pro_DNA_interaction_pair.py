#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob

def main(filename):
    # print(filename)
    dist0_lst = []
    dist3_lst = []
    dist5_lst = []
    angle0_lst = []
    angle3_lst = []
    angle5_lst = []

    with open(filename, 'r') as fin:
        local_dists = []
        icount = 0
        for line in fin:
            words = line.split()
            if len(words) < 2:
                continue
            if line.startswith("CONTACT PAIR"):
                local_dists = [0 for i in range(18)]
                icount = 0
            elif line.startswith("END PAIR"):
                dist0_lst.append([local_dists[0], local_dists[6], local_dists[12]])
                dist5_lst.append([local_dists[2], local_dists[8], local_dists[14]])
                dist3_lst.append([local_dists[4], local_dists[10], local_dists[16]])
                angle0_lst.append([local_dists[1], local_dists[7], local_dists[13]])
                angle5_lst.append([local_dists[3], local_dists[9], local_dists[15]])
                angle3_lst.append([local_dists[5], local_dists[11], local_dists[17]])
            elif words[0] == 'No':
                for i in range(6):
                    local_dists[icount] = -1
                    icount += 1
            else:
                d = float(words[-1])
                local_dists[icount] = d
                icount += 1

    quant_list = [dist0_lst, dist5_lst, dist3_lst, angle0_lst, angle5_lst, angle3_lst]
    quant_names = [r"$Base - C_\alpha$", r"$5' Base - C_\alpha$", r"$3' Base - C_\alpha$",
                   r"$Sugar - Base - C_\alpha$", r"$5' Base - Base - C_\alpha$", r"$3' Base -Base - C_\alpha$"]
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(20, 12.5))
    for i, ql in enumerate(quant_list):
        m, n = i // 3 , i % 3
        X1 = []
        Y1 = []
        X2 = []
        Y2 = []
        for d in ql:
            if d[1] > 0:
                X1.append(d[0])
                Y1.append(d[1])
            if d[2] > 0:
                X2.append(d[0])
                Y2.append(d[2])
        if m == 0:
            datamin, datamax = 4, 20
            databin = 4
        else:
            datamin, datamax = 0, 180
            databin = 45
        x3 = np.arange(datamin, datamax + 1, databin)
        y3 = np.arange(datamin, datamax + 1, databin)
        axes[m, n].plot(X1, Y1, 'r.')
        axes[m, n].plot(X2, Y2, 'g.')
        axes[m, n].plot(x3, y3, 'k--')

        axes[m, n].set_aspect("equal")
        axes[m, n].set_title(quant_names[i])
        axes[m, n].set_xticks(np.arange(datamin, datamax + 1, databin), minor=False)
        axes[m, n].set_yticks(np.arange(datamin, datamax + 1, databin), minor=False)
        axes[m, n].set_xticklabels(np.arange(datamin, datamax + 1, databin), fontsize = 16)
        axes[m, n].set_yticklabels(np.arange(datamin, datamax + 1, databin), fontsize = 16)
        axes[m, n].set_xlim(datamin, datamax)
        axes[m, n].set_ylim(datamin, datamax)
        for axx in ['top', 'bottom', 'left', 'right']:
            axes[m, n].spines[axx].set_linewidth(1.5)
            axes[m, n].xaxis.set_tick_params(length=6, width=1.5)
            axes[m, n].yaxis.set_tick_params(length=6, width=1.5)
    fig.subplots_adjust(hspace=0.2, right=0.79)
    fig_name = filename[:-4] + "_plot.png"
    # fig.savefig(fig_name,  dpi=80)
    plt.show()
 


if __name__ == '__main__':
    import sys
    main(sys.argv[1])

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

def main(filename):
    # print(filename)
    dist0_lst = []
    angle0_lst = []
    angle3_lst = []
    angle5_lst = []

    nullfmt = NullFormatter()

    with open(filename, 'r') as fin:
        local_dists = []
        icount = 0
        bad_flag = 0
        for line in fin:
            words = line.split()
            if len(words) < 2:
                continue
            if line.startswith("CONTACT PAIR"):
                local_dists = [0 for i in range(16)]
                icount = 0
                bad_flag = 0
            elif line.startswith("END PAIR"):
                if bad_flag == 1:
                    continue
                dist0_lst.append([local_dists[0], local_dists[4], local_dists[8], local_dists[12]])
                angle0_lst.append([local_dists[1], local_dists[5], local_dists[9], local_dists[13]])
                angle5_lst.append([local_dists[2], local_dists[6], local_dists[10], local_dists[14]])
                angle3_lst.append([local_dists[3], local_dists[7], local_dists[11], local_dists[15]])
            elif words[0] == 'No':
                bad_flag = 1
            elif words[3] == "NaN":
                bad_flag = 1
            else:
                d = float(words[-1])
                local_dists[icount] = d
                icount += 1

    quant_list = [dist0_lst, angle0_lst, angle5_lst, angle3_lst]
    quant_names = [r"$Base - C_\alpha$", r"$Sugar - Base - C_\alpha$",
                   r"$5' Base - Base - C_\alpha$", r"$3' Base -Base - C_\alpha$"]

    fig = plt.figure(1, figsize=(18, 16))

    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    botomm_h = left_h = left + width + 0.02

    for i, ql in enumerate(quant_list):
        m, n = i // 2 , i % 2

        left_local, width_local = (n + left) * 0.5, 0.5 * width
        bottom_local, height_local = (1 - m + bottom) * 0.5, 0.5 * height
        bottom_h_local, left_h_local = bottom_local + height_local + 0.01, left_local + width_local + 0.01
        rect_scatter = [left_local, bottom_local, width_local, height_local]
        rect_histx = [left_local, bottom_h_local, width_local, 0.1]
        rect_histy = [left_h_local, bottom_local, 0.1, height_local]
        axS_local = plt.axes(rect_scatter)
        axHx_local = plt.axes(rect_histx)
        axHy_local = plt.axes(rect_histy)

        axHx_local.xaxis.set_major_formatter(nullfmt)
        axHy_local.yaxis.set_major_formatter(nullfmt)

        X0 = []
        X1 = []
        X2 = []
        X3 = []
        Y1 = []
        Y2 = []
        Y3 = []
        for d in ql:
            if d[0] < 0:
                continue
            X0.append(d[0])
            if d[1] > 0:
                X1.append(d[0])
                Y1.append(d[1] - d[0])
            if d[2] > 0:
                X2.append(d[0])
                Y2.append(d[2] - d[0])
            if d[3] > 0:
                X3.append(d[0])
                Y3.append(d[3] - d[0])
        if m == 0 and n == 0:
            datamin_x, datamax_x = 4, 16
            databin_x = 4
            datamin_y, datamax_y = -2, 2
            databin_y = 0.5
            lbl_x = r'$d$'
            lbl_y = r'$\Delta d$'
        else:
            datamin_x, datamax_x = 0, 180
            databin_x = 45
            datamin_y, datamax_y = -40, 40
            databin_y = 10
            lbl_x = r'$\theta$'
            lbl_y = r'$\Delta \theta$'
        bins_x = np.arange(datamin_x, datamax_x + (datamax_x - datamin_x)/40, (datamax_x - datamin_x)/40)
        bins_y = np.arange(datamin_y, datamax_y + (datamax_y - datamin_y)/40, (datamax_y - datamin_y)/40)

        axS_local.plot(X1, Y1, 'r,', alpha=0.3)
        axS_local.plot(X2, Y2, 'g,', alpha=0.3)
        axS_local.plot(X3, Y3, 'b,', alpha=0.3)
        axS_local.axhline(y=0, linestyle='--', color='black')

        axHx_local.set_title(quant_names[i])
        axS_local.set_xticks(np.arange(datamin_x, datamax_x + 1, databin_x), minor=False)
        axS_local.set_yticks(np.arange(datamin_y, datamax_y + 1, databin_y), minor=False)
        axS_local.set_xticklabels(np.arange(datamin_x, datamax_x + 1, databin_x), fontsize = 15)
        axS_local.set_yticklabels(np.arange(datamin_y, datamax_y + 1, databin_y), fontsize = 15)
        axS_local.set_xlim(datamin_x, datamax_x)
        axS_local.set_ylim(datamin_y, datamax_y)
        axS_local.set_xlabel(lbl_x, fontsize=20)
        axS_local.set_ylabel(lbl_y, fontsize=20)

        axHx_local.hist(X0, bins=bins_x)
        axHy_local.hist(Y1, bins=bins_y, orientation='horizontal', color='r', alpha=0.3)
        axHy_local.hist(Y2, bins=bins_y, orientation='horizontal', color='g', alpha=0.3)
        axHy_local.axhline(y=0, linestyle='--', color='black')

        axHx_local.set_xlim(axS_local.get_xlim())
        axHy_local.set_ylim(axS_local.get_ylim())
        axHx_local.set_xticks(np.arange(datamin_x, datamax_x + 0.1, databin_x), minor=False)
        axHy_local.set_yticks(np.arange(datamin_y, datamax_y + 0.1, databin_y), minor=False)
        axHx_local.set_yticks(np.arange(0, 16000, 5000), minor=False)
        axHx_local.set_ylabel(r'$(\times 10^3)$', fontsize=16)
        axHy_local.set_xticks(np.arange(0, 61000, 20000), minor=False)
        axHy_local.set_xlabel(r'$(\times 10^3)$', fontsize=16)
        axHx_local.set_yticklabels(['', 5, 10, 15], fontsize = 15)
        axHy_local.set_xticklabels(['', 20, 40, 60], fontsize = 15)
        for axx in ['top', 'bottom', 'left', 'right']:
            axS_local.spines[axx].set_linewidth(1.5)
            axS_local.xaxis.set_tick_params(length=6, width=1.5)
            axS_local.yaxis.set_tick_params(length=6, width=1.5)
            axHx_local.spines[axx].set_linewidth(1.5)
            axHx_local.xaxis.set_tick_params(length=6, width=1.5)
            axHx_local.yaxis.set_tick_params(length=6, width=1.5)
            axHy_local.spines[axx].set_linewidth(1.5)
            axHy_local.xaxis.set_tick_params(length=6, width=1.5)
            axHy_local.yaxis.set_tick_params(length=6, width=1.5)

    fig.subplots_adjust(hspace=0.2, right=0.79)
    fig_name = filename[:-4] + "_mutation_diff_plot.png"
    fig.savefig(fig_name,  dpi=80)
    # plt.show()
 


if __name__ == '__main__':
    import sys
    main(sys.argv[1])

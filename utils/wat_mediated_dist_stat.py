#!/usr/bin/env python

key_names = ['A', 'T', 'G', 'C', 'P', 'S']
huge_repo = {i : [] for i in key_names}

def anaf(fin_name):
    """Read from distance.dat files.
    Keyword Arguments:
    fin_name -- data file name
    """
    with open(fin_name, 'r') as fin:
        for lines in fin:
            line = lines.split()
            if len(line) == 3 and line[0] == "WAT_MED_PRO":
                key0, dist0 = line[1], float(line[2])
                if dist0 < 2.5:
                    continue
            else:
                continue
            huge_repo[key0].append(dist0)

def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab
    import glob
    for filename in glob.glob(r'./*.dat'):
        print(filename)
        anaf(filename)

    for i in key_names:
        out_name = "wat_mediated_pro_"+ i + "_distro.png"
        x = np.array(huge_repo[i])
        num_bins = 100
        # the histogram of data
        n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
        # add a 'best fit' line
        # y = mlab.normpdf(bins, x.mean(), x.std())
        y_sum = sum(n[:])
        y, z = [sum(n[:i])/y_sum for i in range(len(n)+1)], [0.95 for i in range(len(n)+ 1)]
        plt.plot(bins, y, 'r-')
        plt.plot(bins, z, 'r--')
        plt.xlabel(r'distance $\AA$')
        plt.ylabel('Probability')
        plt.xlim(0, 10)
        plt.xticks(np.arange(11))
        plt.title('Histogram of Distance Between Water-Mediated Protein and '+i)
        plt.savefig(out_name, dpi=150)
        plt.clf()

if __name__ == '__main__':
    main()

#!/usr/bin/env python

key_names = ['A', 'T', 'G', 'C', 'P', 'S', 'PRO']
huge_repo = {i : [] for i in key_names}
N = 250

def anaf(fin_name):
    """Read from distance.dat files.
    Keyword Arguments:
    fin_name -- data file name
    """
    global N
    with open(fin_name, 'r') as fin:
        for lines in fin:
            line = lines.split()
            if len(line) == 3 and line[0] == "WAT_PAIR":
                key0, dist0 = line[1], float(line[2])
                if dist0 < 2:
                    print(fin_name, "  ", dist0, "  ", key0)
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

    # for i in key_names:
    #     out_name = "near_wat_"+ i + "_distro.png"
    #     x = np.array(huge_repo[i])
    #     num_bins = 50
    #     # the histogram of data
    #     n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
    #     # add a 'best fit' line
    #     y = mlab.normpdf(bins, x.mean(), x.std())
    #     plt.plot(bins, y, 'r--')
    #     plt.xlabel(r'distance \AA')
    #     plt.ylabel('Probability')
    #     plt.xlim(0, 10)
    #     plt.title('Histogram of Distance from Water to '+i)
    #     plt.savefig(out_name, dpi=150)
    #     plt.clf()

    # for i in key_names:
    #     out_name = "wat_"+ i + "_distro.png"
    #     x = np.array(huge_repo[i])
    #     num_bins = 50
    #     # the histogram of data
    #     n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
    #     # add a 'best fit' line
    #     y = mlab.normpdf(bins, x.mean(), x.std())
    #     plt.plot(bins, y, 'r--')
    #     plt.xlabel(r'distance \AA')
    #     plt.ylabel('Probability')
    #     plt.xlim(0, 40)
    #     plt.title('Histogram of Distance from Water to '+i)
    #     plt.savefig(out_name, dpi=150)
    #     plt.clf()

if __name__ == '__main__':
    main()

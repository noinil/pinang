#!/usr/bin/env python

def b_ds(f):
    """ Calculate length at a given force f (dsDNA):
    Keyword Arguments:
    f -- force
    Worm-Like Chain (WLC)
    """
    Bds = 0.34*10**-9
    kBT = 1.38*10**-23 * 300
    Pds = 48*10**-9
    Sds = 1200*10**-12
    return Bds * (1 - 0.5 * (kBT / (Pds * f * 10**-12))**0.5 + f * 10**-12/Sds )/10**-9

def b_ss(f):
    """ Calculate length at a given force f (ssDNA):
    Keyword Arguments:
    f -- force
    Freely-Jointed Chain (FJC)
    """
    from math import tanh
    Bss = 0.55*10**-9
    kBT = 1.38*10**-23 * 300
    Pss = 0.75*10**-9
    Sss = 720*10**-12
    return Bss * (1 / (tanh(2 * Pss * f * 10**-12/kBT)) - 0.5 * (kBT / (Pss * f * 10**-12)) ) * (1 + f * 10**-12/Sss )/10**-9


def main():
    import numpy as np
    import matplotlib.pyplot as plt

    X = [i*0.1 for i in range(1, 950)]
    Y_ds = [b_ds(i) for i in X]
    Y_ss = [b_ss(i) for i in X]

    plt.xlim(0.01, 0.65)
    plt.ylim(0, 100)
    plt.plot(Y_ds, X)
    plt.plot(Y_ss, X)
    plt.show()

if __name__ == '__main__':
    main()

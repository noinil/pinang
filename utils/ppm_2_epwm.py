#!/usr/bin/env python3

import sys
import numpy as np

def main(ppmname):
    ppm = {}
    with open(ppmname, 'r') as ppm_fin:
        for line in ppm_fin:
            words = line.split()
            if len(words) < 1:
                continue
            if words[0] in ['A:', 'C:', 'G:', 'T:']:
                local_p = []
                for ppm_value in words[1:]:
                    v = float(ppm_value)
                    local_p.append(v)
                ppm[words[0][0]] = local_p[:]

    en_from_p = {}
    en_from_p['A'] = []
    en_from_p['C'] = []
    en_from_p['G'] = []
    en_from_p['T'] = []
    len_DNA = len(ppm['A'])
    for i in range(len_DNA):
        pA, pC, pG, pT = ppm['A'][i], ppm['C'][i], ppm['G'][i], ppm['T'][i]
        eA, eC, eG, eT = -np.log(pA/0.25), -np.log(pC/0.25), -np.log(pG/0.25), -np.log(pT/0.25)
        e_ave = (eA + eC + eG + eT) / 4
        eA -= e_ave
        eC -= e_ave
        eG -= e_ave
        eT -= e_ave
        en_from_p['A'].append(eA)
        en_from_p['C'].append(eC)
        en_from_p['G'].append(eG)
        en_from_p['T'].append(eT)

    # Output
    outstr = " {0:8.4f} "
    outfilename = ppmname[:-3] + "pwm"
    outfile = open(outfilename, 'w')
    for bname in ['A', 'C', 'G', 'T']:
        outfile.write(bname + ": ")
        for v in en_from_p[bname]:
            outfile.write(outstr.format(v))
        outfile.write(' \n')
    outfile.close()


if __name__ == '__main__':
    ppm_fname = sys.argv[1]
    main(ppm_fname)

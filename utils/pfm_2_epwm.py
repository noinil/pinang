#!/usr/bin/env python3

import sys
import numpy as np

def main(pfmname):
    pfm = {}
    with open(pfmname, 'r') as pwm_fin:
        for line in pwm_fin:
            words = line.split()
            if len(words) < 1:
                continue
            if words[0] in ['A:', 'C:', 'G:', 'T:']:
                local_p = []
                for pfm_value in words[1:]:
                    v = float(pfm_value)
                    local_p.append(v)
                pfm[words[0][0]] = local_p[:]

    en_from_p = {}
    en_from_p['A'] = []
    en_from_p['C'] = []
    en_from_p['G'] = []
    en_from_p['T'] = []
    len_DNA = len(pfm['A'])
    for i in range(len_DNA):
        pA, pC, pG, pT = pfm['A'][i], pfm['C'][i], pfm['G'][i], pfm['T'][i]
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
    outfilename = pfmname[:-3] + "pwm"
    outfile = open(outfilename, 'w')
    for bname in ['A', 'C', 'G', 'T']:
        outfile.write(bname + ": ")
        for v in en_from_p[bname]:
            outfile.write(outstr.format(v))
        outfile.write(' \n')
    outfile.close()


if __name__ == '__main__':
    pfm_fname = sys.argv[1]
    main(pfm_fname)

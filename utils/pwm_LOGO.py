#!/usr/bin/env python3

import sys
import numpy as np
import operator

def main(pwm_file_name):
    # ---------- read in probability matrix ----------
    p_pwm = {}
    with open(pwm_file_name, 'r') as pwm_fin:
        for line in pwm_fin:
            words = line.split()
            if len(words) < 1:
                continue
            if words[0] in ['A:', 'C:', 'G:', 'T:']:
                local_p = []
                for pwm_value in words[1:]:
                    v = float(pwm_value)
                    local_p.append(v)
                p_pwm[words[0][0]] = local_p[:]

    len_A = len(p_pwm['A'])
    len_C = len(p_pwm['C'])
    len_G = len(p_pwm['G'])
    len_T = len(p_pwm['T'])
    if len_A != len_C or len_C != len_G or len_G != len_T:
        print("Inconsistent length of PWM!")
        return
    len_DNA = len_A

    # ---------- Plotting ----------
    svg_header_line = '''<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {0} {1}" version="1.1">
'''
    svg_title_line = '<title>{0} LOGO</title> \n'
    svg_desc_line = '<desc>LOGO of motif: {0}, generated from PPM</desc> \n'
    svg_pwm_text_char_line = '<text x="{0}" y="{1:4.2f}" fill="{2}" style="font-weight: bold; font-family: \'{3}\'; font-size: {4:4.2f}pt" textLength="20" lengthAdjust="spacingAndGlyphs">{5}</text> \n'
    svg_pwm_height = 40
    svg_pwm_width = 20
    svg_color_dict = {'A': "#4cdef5", 'C': "#a4d555", 'G': "#ff5992", 'T': "#841983"}

    svg_file_name = pwm_file_name[:-3] + 'svg'
    svg_file = open(svg_file_name, 'w')
    svg_file.write(svg_header_line.format((len_DNA + 2) * svg_pwm_width, 2 * svg_pwm_height))
    svg_file.write(svg_title_line.format(pwm_file_name[:-4]))
    svg_file.write(svg_desc_line.format(pwm_file_name[:-4]))

    p_base_dict_local = {}
    x, y = 20, 20
    for i in range(len_DNA):
        y = 0
        pA, pC, pG, pT = p_pwm['A'][i], p_pwm['C'][i], p_pwm['G'][i], p_pwm['T'][i]
        eA, eC, eG, eT = np.log2(pA/0.25), np.log2(pC/0.25), np.log2(pG/0.25), np.log2(pT/0.25)
        IC = pA * eA + pC * eC + pG * eG + pT * eT
        p_base_dict_local['A'] = pA
        p_base_dict_local['C'] = pC
        p_base_dict_local['G'] = pG
        p_base_dict_local['T'] = pT
        sorted_p_base = sorted(p_base_dict_local.items(), key=operator.itemgetter(1), reverse=True)
        # print(sorted_p_base)
        y += svg_pwm_height * (1 - IC / 2.0)
        for j in sorted_p_base:
            base_name, p = j[0], j[1]
            dy = p * IC * svg_pwm_height / 2
            y += dy
            svg_file.write(svg_pwm_text_char_line.format(x, y, svg_color_dict[base_name], "Arial", dy, base_name))
        svg_file.write(" \n")
        x += svg_pwm_width

    svg_file.write('</svg>\n')
    svg_file.close()

if __name__ == '__main__':
    try:
        pwm_file_name = sys.argv[1]
    except:
        print("Usage: ", sys.argv[0], " xxx.ppm")
        exit()
    main(pwm_file_name)

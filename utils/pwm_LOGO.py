#!/usr/bin/env python3

import sys
import numpy as np
import operator

def main(pwm_file_name, reverse_option):
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

    if reverse_option:
        p_new_pwm_a = list(reversed(p_pwm['T']))
        p_new_pwm_t = list(reversed(p_pwm['A']))
        p_new_pwm_g = list(reversed(p_pwm['C']))
        p_new_pwm_c = list(reversed(p_pwm['G']))
        p_pwm['A'] = p_new_pwm_a
        p_pwm['C'] = p_new_pwm_c
        p_pwm['G'] = p_new_pwm_g
        p_pwm['T'] = p_new_pwm_t

    # ---------- Plotting ----------
    svg_header_line = '''<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {0} {1}" version="1.1">
'''
    svg_title_line = '<title>{0} LOGO</title> \n'
    svg_desc_line = '<desc>LOGO of motif: {0}, generated with PINANG.  Copyright: Cheng Tan, https://c-tan.com </desc> \n'
    svg_pwm_text_char_line = '<text x="{0}" y="{1:4.2f}" fill="{2}" style="font-weight: bold; font-family: \'{3}\'; font-size: {4:4.2f}pt" textLength="20" lengthAdjust="spacingAndGlyphs">{5}</text> \n'
    svg_pwm_straight_vline = '<line x1="20" y1="0" x2="20" y2="40" style="stroke:rgb(0,0,0);stroke-width:0.5" /> \n'
    svg_pwm_straight_vtic2 = '<line x1="17" y1="0" x2="20" y2="0" style="stroke:rgb(0,0,0);stroke-width:0.5" /> \n'
    svg_pwm_straight_vtic1 = '<line x1="17" y1="20" x2="20" y2="20" style="stroke:rgb(0,0,0);stroke-width:0.5" /> \n'
    svg_pwm_straight_vtic0 = '<line x1="17" y1="40" x2="20" y2="40" style="stroke:rgb(0,0,0);stroke-width:0.5" /> \n'
    svg_pwm_straight_hline = '<line x1="20" y1="40" x2="{0}" y2="40" style="stroke:rgb(0,0,0);stroke-width:0.5" /> \n'
    svg_pwm_height = 40
    svg_pwm_width = 20
    svg_color_dict = {'A': "#b50804", 'C': "#83bcc3", 'G': "#eaa612", 'T': "#3e4079"}

    svg_file_name = pwm_file_name[:-3] + 'svg'
    svg_file = open(svg_file_name, 'w')
    svg_file.write(svg_header_line.format((len_DNA + 2) * svg_pwm_width, 2 * svg_pwm_height))
    svg_file.write(svg_title_line.format(pwm_file_name[:-4]))
    svg_file.write(svg_desc_line.format(pwm_file_name[:-4]))
    svg_file.write(svg_pwm_straight_vline)
    svg_file.write(svg_pwm_straight_vtic2)
    svg_file.write(svg_pwm_straight_vtic1)
    svg_file.write(svg_pwm_straight_vtic0)
    svg_file.write(svg_pwm_straight_hline.format((len_DNA + 1) * svg_pwm_width))

    p_base_dict_local = {}
    x, y = 20, 20
    for i in range(len_DNA):
        y = 0
        fA, fC, fG, fT = p_pwm['A'][i], p_pwm['C'][i], p_pwm['G'][i], p_pwm['T'][i]
        f_sum = fA + fC + fG + fT
        pA, pC, pG, pT = fA / f_sum, fC / f_sum, fG / f_sum, fT / f_sum, 
        eA, eC, eG, eT = np.log2(pA/0.25), np.log2(pC/0.25), np.log2(pG/0.25), np.log2(pT/0.25)
        IC = pA * eA + pC * eC + pG * eG + pT * eT
        p_base_dict_local['A'] = pA
        p_base_dict_local['C'] = pC
        p_base_dict_local['G'] = pG
        p_base_dict_local['T'] = pT
        sorted_p_base = sorted(p_base_dict_local.items(), key=operator.itemgetter(1), reverse=True)
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
    import argparse
    def parse_arguments():
        parser = argparse.ArgumentParser(description='Plot energy LOGO from PWM/PFM file.')
        parser.add_argument('-r', '--reverse', action="store_true", help="Reverse sequence PWM")
        parser.add_argument('pwmFileName', type=str, help="Position Probability/Frequency Matrix")
        return parser.parse_args()
    args = parse_arguments()
    main(args.pwmFileName, args.reverse)

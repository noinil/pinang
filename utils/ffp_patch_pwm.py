#!/usr/bin/env python3

import sys

def main(ffp_file_name, pwm_file_name):

    # -------------------- pwm read in --------------------
    i_pair = {}                 # convert resid ID to i of pairs
    ene_pair = {}               # convert i of pairs to energies
    ip_ene = 0
    ip_count = {}
    strandA, strandB = [], []
    with open(pwm_file_name, 'r') as pwmf:
        for line in pwmf:
            words = line.split()
            if len(words) < 2:
                continue
            a, b = int(words[0]), int(words[1])
            eA, eC, eG, eT = float(words[2]), float(words[3]), float(words[4]), float(words[5])
            i_pair[a] = ip_ene
            i_pair[b] = ip_ene
            strandA.append(a)
            strandB.append(b)
            ene_pair[ip_ene] = (eA, eC, eG, eT)
            ip_count[ip_ene] = 0
            ip_ene += 1

    # -------------------- FFP rewrite --------------------
    new_ffp_name = ffp_file_name[:-4] + "_patched.ffp"
    new_ffpf = open(new_ffp_name, 'w')

    new_ffpf_section = "[ protein-DNA seq-specific w/ PWM ]  {0:6d} \n"
    new_ffpf_desc = "#  pro_i         r_0 angle_NC  angle_0 angle_53    sigma      phi     ene_A     ene_C     ene_G     ene_T \n"
    new_ffpf_pwm = "{0:8d} {1:11.6f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f} {7:9.6f} {8:9.6f} {9:9.6f} {10:9.6f} \n"
    rewrite_flag = 0
    interaction_pair_num = 0
    ffp_data = []
    with open(ffp_file_name, 'r') as ffpf:
        for line in ffpf:
            words = line.split()
            if line.startswith('[ protein-DNA seq-specific ]'):
                interaction_pair_num = int(words[-1])  
                rewrite_flag = interaction_pair_num 
                continue
            if rewrite_flag == 0:
                new_ffpf.write(line)
            elif not line.startswith('#'):
                i_pro, i_dna = int(words[0]), int(words[1])
                r0, a_nc, a_0, a_53 = float(words[2]), float(words[3]), float(words[4]), float(words[5])
                sigma, phi = float(words[6]), float(words[7])
                if i_dna in i_pair.keys():
                    ip_count[i_pair[i_dna]] += 1
                ffp_data.append((i_pro, i_dna, r0, a_nc, a_0, a_53, sigma, phi))
                rewrite_flag -= 1


    for k, v in ene_pair.items():
        i = ip_count[k]
        if i > 0:
            eA, eC, eG, eT = v[0]/i, v[1]/i, v[2]/i, v[3]/i
        ene_pair[k] = (eA, eC, eG, eT)
        
    new_ffpf.write(new_ffpf_section.format(interaction_pair_num))
    new_ffpf.write(new_ffpf_desc)
    for d in ffp_data:
        i_pro, i_dna = d[0], d[1]
        r0, a_nc, a_0, a_53 = d[2], d[3], d[4], d[5]
        sigma, phi = d[6], d[7]
        v = ene_pair[i_pair[i_dna]]
        if i_dna in strandA:
            eA, eC, eG, eT = v[0], v[1], v[2], v[3]
        elif i_dna in strandB:
            eA, eC, eG, eT = v[3], v[2], v[1], v[0]
        else:
            eA, eC, eG, eT = 0, 0, 0, 0
        new_ffpf.write(new_ffpf_pwm.format(i_pro, r0, a_nc, a_0, a_53, sigma, phi, eA, eC, eG, eT))

    new_ffpf.close()


if __name__ == '__main__':
    try:
        ffp_file_name = sys.argv[1]
        pwm_file_name = sys.argv[2]
    except:
        print(" Usage: ", sys.argv[0], " xxx.ffp xxx.pwm")
    main(ffp_file_name, pwm_file_name)

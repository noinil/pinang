#!/usr/bin/env python3

import numpy as np

def read_pfm_to_pwm(filename):
    pfm = {}
    with open(filename, 'r') as pfm_fin:
        for line in pfm_fin:
            words = line.split()
            if len(words) < 1:
                continue
            w0 = words[0]
            if w0 in ['A:', 'C:', 'G:', 'T:']:
                local_list = []
                for pfm_value in words[1:]:
                    v = float(pfm_value)
                    local_list.append(v)
                pfm[w0] = local_list[:]
            elif w0 in {"CHAIN_A:", "CHAIN_B:"}:
                local_list = []
                for dna_res_id in words[1:]:
                    k = int(dna_res_id)
                    local_list.append(k)
                pfm[w0] = local_list[:]
    ffp_pfm_T = np.array([pfm['A:'], pfm['C:'], pfm['G:'], pfm['T:']])
    ffp_pfm_col_sum = np.sum(ffp_pfm_T, axis=0)
    ffp_ppm_T = ffp_pfm_T / ffp_pfm_col_sum
    ffp_pwm_T_0 = -np.log(ffp_ppm_T)
    ffp_pwm_col_sum = np.sum(ffp_pwm_T_0, axis=0)
    ffp_pwm_T = ffp_pwm_T_0 - ffp_pwm_col_sum / 4
    return (ffp_pwm_T.T, pfm['CHAIN_A:'], pfm['CHAIN_B:'])

def read_ppm_to_pwm(filename):
    ppm = {}
    with open(filename, 'r') as ppm_fin:
        for line in ppm_fin:
            words = line.split()
            if len(words) < 1:
                continue
            w0 = words[0]
            if w0 in ['A:', 'C:', 'G:', 'T:']:
                local_list = []
                for ppm_value in words[1:]:
                    v = float(ppm_value)
                    local_list.append(v)
                ppm[w0] = local_list[:]
            elif w0 in {"CHAIN_A:", "CHAIN_B:"}:
                local_list = []
                for dna_res_id in words[1:]:
                    k = int(dna_res_id)
                    local_list.append(k)
                ppm[w0] = local_list[:]
    ffp_ppm_T = np.array([ppm['A:'], ppm['C:'], ppm['G:'], ppm['T:']])
    ffp_pwm_T_0 = -np.log(ffp_ppm_T)
    ffp_pwm_col_sum = np.sum(ffp_pwm_T_0, axis=0)
    ffp_pwm_T = ffp_pwm_T_0 - ffp_pwm_col_sum / 4
    return (ffp_pwm_T.T, ppm['CHAIN_A:'], ppm['CHAIN_B:'])


def main(pwm_option, pwm_name, ffp_name):
    if pwm_option == '-p':
        ffp_pwm, chain_A_id, chain_B_id = read_ppm_to_pwm(pwm_name)
    elif pwm_option == '-f':
        ffp_pwm, chain_A_id, chain_B_id = read_pfm_to_pwm(pwm_name)
    pwm_len = len(chain_A_id)

    # rough checking ......
    def chain_id_check(a):
        for i in range(1, len(a)):
            if a[i] - a[i - 1] != 3:
                print(" Chain ID error!  Please check!")
                exit()
    chain_id_check(chain_A_id)
    chain_id_check(chain_B_id[::-1])
    if len(chain_B_id) != pwm_len:
        print(" Inconsistent length of CHAIN ids!  Please check!")
        exit()
    # ------------------------------------------------------------

    # -------------------- FFP rewrite --------------------
    ip_count = np.zeros((pwm_len, 1), dtype=np.int)  # count for decomposition

    new_ffp_name = ffp_name[:-4] + "_patched.ffp"
    new_ffp = open(new_ffp_name, 'w')

    new_ffpf_section = "[ protein-DNA seq-specific w/ PWM ]  {0:6d} \n"
    new_ffpf_desc = "#  pro_i         r_0 angle_NC  angle_0 angle_53     ene_A     ene_C     ene_G     ene_T    sigma      phi \n"
    new_ffpf_pwm = "{0:8d} {1:11.6f} {2:8.3f} {3:8.3f} {4:8.3f} {7:9.6f} {8:9.6f} {9:9.6f} {10:9.6f} {5:8.3f} {6:8.3f} \n"
    rewrite_flag = 0
    interaction_pair_num = 0
    ffp_data = []
    with open(ffp_name, 'r') as ffp:
        for line in ffp:
            words = line.split()
            if line.startswith('[ protein-DNA seq-specific ]'):
                interaction_pair_num = int(words[-1])  
                rewrite_flag = interaction_pair_num 
                continue
            if rewrite_flag == 0:
                new_ffp.write(line)
            elif not line.startswith('#'):
                i_pro, i_dna = int(words[0]), int(words[1])
                r0, a_nc, a_0, a_53, sigma, phi = float(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6]), float(words[7])
                try:
                    if i_dna in chain_A_id:
                        k = chain_A_id.index(i_dna)
                        j = 1
                    else:
                        k = chain_B_id.index(i_dna)
                        j = -1
                    ip_count[k] += 1
                except:
                    print(" Index error in CHAIN_id s!   Please Check!")
                    exit()
                ffp_data.append((i_pro, i_dna, r0, a_nc, a_0, a_53, sigma, phi, k, j))
                rewrite_flag -= 1

    for i, k in enumerate(ip_count):
        if k == 0:
            ip_count[i] = -1
    ffp_pwm_decomposed = ffp_pwm / ip_count
        
    new_ffp.write(new_ffpf_section.format(interaction_pair_num))
    new_ffp.write(new_ffpf_desc)
    for d in ffp_data:
        i_pro, i_dna = d[0], d[1]
        r0, a_nc, a_0, a_53, sigma, phi = d[2], d[3], d[4], d[5], d[6], d[7]
        k, j = d[8], d[9]
        v = ffp_pwm_decomposed[k]
        eA, eC, eG, eT = tuple(v[::j])
        new_ffp.write(new_ffpf_pwm.format(i_pro, r0, a_nc, a_0, a_53, sigma, phi, eA, eC, eG, eT))
    new_ffp.close()


if __name__ == '__main__':
    import sys
    pwm_option = sys.argv[1]
    if pwm_option not in {'-p', '-f'}:
        print("Usage: ", sys.argv[0], " [-p pro.ppm] [-f pro.pfm] -P pro.ffp")
        exit()
    pwm_name = sys.argv[2]
    ffp_name = sys.argv[4]
    main(pwm_option, pwm_name, ffp_name)

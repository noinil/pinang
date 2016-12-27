#!/usr/bin/env python3

def main():
    en_pwm = {}
    with open('en.pwm', 'r') as pwm_fin:
        for line in pwm_fin:
            words = line.split()
            if len(words) < 1:
                continue
            if words[0] in ['A:', 'C:', 'G:', 'T:']:
                local_list = []
                for pwm_value in words[1:]:
                    v = float(pwm_value)
                    local_list.append(v)
                en_pwm[words[0][0]] = local_list[:]
            elif words[0] == "CHAIN_A:":
                local_list = []
                for dna_res_id in words[1:]:
                    k = int(dna_res_id)
                    # local_list.append(k * 3 - 1)
                    local_list.append(k)
                en_pwm["chain_A"] = local_list[:]
            elif words[0] == "CHAIN_B:":
                local_list = []
                for dna_res_id in words[1:]:
                    k = int(dna_res_id)
                    # local_list.append(k * 3 - 2)
                    local_list.append(k)
                en_pwm["chain_B"] = local_list[:]

    out_str = " {0:6d} {1:6d} {2:12.8f} {3:12.8f} {4:12.8f} {5:12.8f} "
    len_DNA = len(en_pwm['A'])
    for i in range(len_DNA):
        print(out_str.format(en_pwm['chain_A'][i], en_pwm['chain_B'][i], en_pwm['A'][i], en_pwm['C'][i], en_pwm['G'][i], en_pwm['T'][i]))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

def main(ifname):
    pro_iunit_junit_list = []
    all_lines = []
    line_output_index = []
    with open(ifname, 'r') as fin:
        for line in fin:
            words = line.split()
            if len(words) < 1:
                all_lines.append(line)
                line_output_index.append(0)
            elif line.startswith("<<<<") or line.startswith("**") or line.startswith('>>>>'):
                all_lines.append(line)
                line_output_index.append(0)
                # pass
            elif words[0] in ['bond', 'angl', 'aicg13', 'dihd', 'aicgdih', 'aicg14']:
                if "pp" in line:
                    all_lines.append(line)
                    iunit, junit = int(words[2]), int(words[3])
                    ij_pair = (iunit, junit)
                    line_output_index.append( ij_pair )
                    if ij_pair not in pro_iunit_junit_list:
                        pro_iunit_junit_list.append(ij_pair)
            elif words[0] == 'contact':
                f = float(words[-2])
                if f < 0.0001:
                    continue
                iunit, junit = int(words[2]), int(words[3])
                ij_pair = (iunit, junit)
                line_output_index.append( ij_pair )
                all_lines.append(line)
                if ij_pair not in pro_iunit_junit_list:
                    pro_iunit_junit_list.append(ij_pair)

    ninfo_out_name_list = ["pro_{0[0]:0>2d}_{0[1]:0>2d}.ninfo".format(ij) for ij in pro_iunit_junit_list]
    ninfo_out_files = [open(i, 'w') for i in ninfo_out_name_list]
    for i, line in enumerate(all_lines):
        out_status = line_output_index[i]
        if out_status == 0:
            for out_f in ninfo_out_files:
                out_f.write(line)
        else:
            out_index = pro_iunit_junit_list.index(out_status)
            out_f = ninfo_out_files[out_index]
            out_f.write(line)

    for f in ninfo_out_files:
        f.close()

if __name__ == '__main__':
    import sys
    ninfo_file_name = sys.argv[1]
    main(ninfo_file_name)

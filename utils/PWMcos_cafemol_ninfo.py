#!/usr/bin/env python3

import numpy as np
import MDAnalysis 

def main(PDB_name, PFM_name, flag_psf_output, flag_native_contact_list, gamma=1.0, epsilon_prime=0.0, cutoff=4.0):
    # ==========================
    # Variables can be modified:
    # ==========================
    chainID_DNA_a = 'A'
    chainID_DNA_b = 'B'
    protein_name = PFM_name[:-4]

    PRO_DNA_AA_CUTOFF = cutoff

    # ============================================
    # Constants for CG particle masses and charges
    # ============================================
    std_base_mass = {'A': 134.1, 'G': 150.1, 'C': 110.1, 'T': 125.1}
    std_prot_mass = {'ALA' :  71.09,
                     'ARG' :  156.19,
                     'ASN' :  114.11,
                     'ASP' :  115.09,
                     'CYS' :  103.15,
                     'GLN' :  128.14,
                     'GLU' :  129.12,
                     'GLY' :  57.05,
                     'HIS' :  137.14,
                     'ILE' :  113.16,
                     'LEU' :  113.16,
                     'LYS' :  128.17,
                     'MET' :  131.19,
                     'PHE' :  147.18,
                     'PRO' :  97.12,
                     'SER' :  87.08,
                     'THR' :  101.11,
                     'TRP' :  186.21,
                     'TYR' :  163.18,
                     'VAL' :  99.14,
                     'OTH' :  100.0 }
    std_prot_charge = {'ALA' : 0.0,
                       'ARG' : 1.0,
                       'ASN' : 0.0,
                       'ASP' : -1.0,
                       'CYS' : 0.0,
                       'GLN' : 0.0,
                       'GLU' : -1.0,
                       'GLY' : 0.0,
                       'HIS' : 0.0,
                       'ILE' : 0.0,
                       'LEU' : 0.0,
                       'LYS' : 1.0,
                       'MET' : 0.0,
                       'PHE' : 0.0,
                       'PRO' : 0.0,
                       'SER' : 0.0,
                       'THR' : 0.0,
                       'TRP' : 0.0,
                       'TYR' : 0.0,
                       'VAL' : 0.0,
                       'OTH' : 0.0}

    # ============================================
    # Read in structural informatino from PDB file
    # ============================================
    u = MDAnalysis.Universe(PDB_name)

    selstr_DNA_a = "nucleic and segid {0}".format(chainID_DNA_a)
    selstr_DNA_b = "nucleic and segid {0}".format(chainID_DNA_b)
    selstr_protein = "protein and not name H*"

    selstr_P = "(resid {0} and (name P or name OP* or name O5')) or (resid {1} and name O3')"
    selstr_S = "resid {0} and (name C5' or name C4' or name C3' or name C2' or name C1' or name O4' or name O2')"
    selstr_B = "resid {0} and not (name C5' or name C4' or name C3' or name C2' or name C1' or name O4' or name O2' or name O3' or name O5' or name OP* or name P or name H*)"

    sel_dna_a = u.select_atoms(selstr_DNA_a)
    sel_dna_b = u.select_atoms(selstr_DNA_b)
    sel_protein = u.select_atoms(selstr_protein)

    def cg_dna_top(dna_atom_group):
        """Translate atomistic DNA structure into coarse-grained topology.
        Input: dna_atom_group - Atomistic DNA atom/residue group.
        Output: Coarse-grained DNA particles.
        """
        cg_dna_particle_num = len(dna_atom_group.residues) * 3 - 1
        cg_dna_coors = np.empty([cg_dna_particle_num, 3])
        cg_dna_particle_resname = []
        cg_dna_particle_name = []
        cg_dna_particle_charge = []
        cg_dna_particle_mass = []
        resid_list = list(dna_atom_group.residues.resids)
        for i, j in enumerate(resid_list):
            tmp_resname = dna_atom_group.residues[i].resname
            # Phosphate
            if i > 0:
                res_P = dna_atom_group.select_atoms(selstr_P.format(j, j-1))
                cg_dna_coors[i * 3 - 1] = res_P.center_of_mass()
                cg_dna_particle_name.append('DP')
                cg_dna_particle_charge.append(-0.6)
                cg_dna_particle_mass.append(94.97)
                cg_dna_particle_resname.append(tmp_resname)
            # Sugar
            res_S = dna_atom_group.select_atoms(selstr_S.format(j))
            cg_dna_coors[i * 3] = res_S.center_of_mass()
            cg_dna_particle_name.append('DS')
            cg_dna_particle_charge.append(0.0)
            cg_dna_particle_mass.append(83.11)
            cg_dna_particle_resname.append(tmp_resname)
            # Base
            res_B = dna_atom_group.select_atoms(selstr_B.format(j))
            cg_dna_coors[i * 3 + 1] = res_B.center_of_mass()
            cg_dna_particle_name.append('DB')
            cg_dna_particle_charge.append(0.0)
            cg_dna_particle_mass.append(std_base_mass[tmp_resname[-1]])
            cg_dna_particle_resname.append(tmp_resname)
        return (cg_dna_coors, cg_dna_particle_name, cg_dna_particle_charge, cg_dna_particle_mass, cg_dna_particle_resname)

    def cg_pro_top(pro_atom_group):
        """Translate atomistic protein structure into coarse-grained topology.
        Input: pro_atom_group - Atomistic protein atom/residue group.
        Output: Coarse-grained C_alpha particles.
        """
        sel_tmp_all_ca = pro_atom_group.select_atoms("name CA")
        cg_pro_coors = sel_tmp_all_ca.positions
        cg_pro_particle_name = list(sel_tmp_all_ca.names)
        cg_pro_particle_resname = list(sel_tmp_all_ca.resnames)
        cg_pro_particle_segid = list(sel_tmp_all_ca.segids)
        cg_pro_particle_charge = []
        cg_pro_particle_mass = []
        for i, tmp_resname in enumerate(cg_pro_particle_resname):
            cg_pro_particle_charge.append(std_prot_charge[tmp_resname])
            cg_pro_particle_mass.append(std_prot_mass[tmp_resname])
        return (cg_pro_coors, cg_pro_particle_name, cg_pro_particle_segid, cg_pro_particle_charge, cg_pro_particle_mass, cg_pro_particle_resname)

    cg_dna_a_coors, cg_dna_a_p_name, cg_dna_a_p_charge, cg_dna_a_p_mass, cg_dna_a_p_resname = cg_dna_top(sel_dna_a)
    cg_dna_b_coors, cg_dna_b_p_name, cg_dna_b_p_charge, cg_dna_b_p_mass, cg_dna_b_p_resname = cg_dna_top(sel_dna_b)
    cg_pro_coors, cg_pro_p_name, cg_pro_p_segid, cg_pro_p_charge, cg_pro_p_mass, cg_pro_p_resname = cg_pro_top(sel_protein)

    # Assign CG particle ID and residue ID
    cg_dna_a_p_ID = []
    cg_dna_b_p_ID = []
    cg_pro_p_ID = []

    cg_dna_a_r_ID = []
    cg_dna_b_r_ID = []
    cg_pro_r_ID = []

    cg_dna_a_p_num = cg_dna_a_coors.shape[0]
    cg_dna_b_p_num = cg_dna_b_coors.shape[0]
    cg_pro_p_num = cg_pro_coors.shape[0]
    cg_top_total_num = cg_dna_a_p_num + cg_dna_b_p_num + cg_pro_p_num

    cg_dna_a_r_num = len(sel_dna_a.residues)
    cg_dna_b_r_num = len(sel_dna_b.residues)
    cg_pro_r_num = len(sel_protein.residues)

    for i in range(cg_dna_a_p_num):
        cg_dna_a_p_ID.append(i + 1)
        cg_dna_a_r_ID.append((i + 1) // 3 + 1)
    for i in range(cg_dna_b_p_num):
        cg_dna_b_p_ID.append(i + 1 + cg_dna_a_p_num)
        cg_dna_b_r_ID.append((i + 1) // 3 + 1 + cg_dna_a_r_num)
    for i in range(cg_pro_p_num):
        cg_pro_p_ID.append(i + 1 + cg_dna_a_p_num + cg_dna_b_p_num)
        cg_pro_r_ID.append(i + 1 + cg_dna_a_r_num + cg_dna_b_r_num)


    def output_psf(protein_name):
        """Output psf file for protein-DNA complex.
        """
        PSF_HEAD_STR = "PSF CMAP \n\n"
        PSF_TITLE_STR0 = "      3 !NTITLE \n"
        PSF_TITLE_STR1 = "REMARKS PSF file created with PWMcos tools. \n"
        PSF_TITLE_STR2 = "REMARKS Protein: {0}  \n"
        PSF_TITLE_STR5 = "REMARKS ======================================== \n"
        PSF_TITLE_STR6 = "       \n"
        PSF_TITLE_STR = PSF_TITLE_STR0 + PSF_TITLE_STR1 + PSF_TITLE_STR2 + PSF_TITLE_STR5 + PSF_TITLE_STR6
        PSF_ATOM_TITLE = " {atom_num:>6d} !NATOM \n"
        PSF_ATOM_LINE = " {atom_ser:>6d} {seg_id:>3} {res_ser:>5d} {res_name:>3} {atom_name:>3} {atom_type:>5}  {charge:>10.6f}  {mass:>10.6f}          0 \n"

        psf_file_name = protein_name + "_dnp_cg.psf"
        psf_file = open(psf_file_name, 'w')
        psf_file.write(PSF_HEAD_STR)
        psf_file.write(PSF_TITLE_STR.format(protein_name))
        psf_file.write(PSF_ATOM_TITLE.format(atom_num = cg_top_total_num))

        for i in range(cg_dna_a_p_num):
            psf_file.write(PSF_ATOM_LINE.format(atom_ser = i + 1,
                                                seg_id = 'a',
                                                res_ser = (i + 1) // 3 + 1,
                                                res_name = cg_dna_a_p_resname[i],
                                                atom_name = cg_dna_a_p_name[i],
                                                atom_type = cg_dna_a_p_name[i][-1],
                                                charge = cg_dna_a_p_charge[i],
                                                mass = cg_dna_a_p_mass[i]))
        for i in range(cg_dna_b_p_num):
            psf_file.write(PSF_ATOM_LINE.format(atom_ser = i + 1 + cg_dna_a_p_num,
                                                seg_id = 'b',
                                                res_ser = (i + 1) // 3 + 1 + cg_dna_a_r_num,
                                                res_name = cg_dna_b_p_resname[i],
                                                atom_name = cg_dna_b_p_name[i],
                                                atom_type = cg_dna_b_p_name[i][-1],
                                                charge = cg_dna_b_p_charge[i],
                                                mass = cg_dna_b_p_mass[i]))
        for i in range(cg_pro_p_num):
            psf_file.write(PSF_ATOM_LINE.format(atom_ser = i + 1 + cg_dna_a_p_num + cg_dna_b_p_num,
                                                seg_id = cg_pro_p_segid[i],
                                                res_ser = i + 1 + cg_dna_a_r_num + cg_dna_b_r_num,
                                                res_name = cg_pro_p_resname[i],
                                                atom_name = cg_pro_p_name[i],
                                                atom_type = 'C',
                                                charge = cg_pro_p_charge[i],
                                                mass = cg_pro_p_mass[i]))
        psf_file.close()

    if flag_psf_output:
        output_psf(protein_name)

    # ===================================
    # Find out native contacts and output
    # ===================================

    def aagroup_dna(dna_atom_group):
        """Divide atomistic DNA structure into PSB groups.
        Input: dna_atom_group - Atomistic DNA atom/residue group.
        Output: DNA atom groups: P, S, B.
        """
        group_dna_coors = []
        resid_list = list(dna_atom_group.residues.resids)
        for i, j in enumerate(resid_list):
            # Phosphate
            if i > 0:
                res_P = dna_atom_group.select_atoms(selstr_P.format(j, j-1))
                coors_res_P = list(res_P.positions)
                group_dna_coors.append(coors_res_P)
            # Sugar
            res_S = dna_atom_group.select_atoms(selstr_S.format(j))
            coors_res_S = list(res_S.positions)
            group_dna_coors.append(coors_res_S)
            # Base
            res_B = dna_atom_group.select_atoms(selstr_B.format(j))
            coors_res_B = list(res_B.positions)
            group_dna_coors.append(coors_res_B)
        return group_dna_coors

    def aagroup_pro(pro_atom_group):
        """Divide atomistic protein structure into aa groups.
        Input: pro_atom_group - Atomistic protein atom/residue group.
        Output: Amino acid groups.
        """
        group_pro_coors = []
        for aa_res in pro_atom_group.residues:
            coors_aa_res = list(aa_res.atoms.positions)
            group_pro_coors.append(coors_aa_res)
        return group_pro_coors

    aagroup_dna_a_coors = aagroup_dna(sel_dna_a)
    aagroup_dna_b_coors = aagroup_dna(sel_dna_b)
    aagroup_pro_coors = aagroup_pro(sel_protein)

    def min_distance_residues(coors_r1, coors_r2):
        mind = 100000
        if len(coors_r1) * len(coors_r2) == 0:
            return
        for c1 in coors_r1:
            for c2 in coors_r2:
                d = np.linalg.norm(c1 - c2)
                if d < mind:
                    mind = d
        return mind
    def vec_angle(vec1, vec2):
        n_v1 = np.linalg.norm(vec1)
        n_v2 = np.linalg.norm(vec2)
        return np.arccos(np.clip(np.dot(vec1, vec2) / n_v1 / n_v2, -1.0, 1.0)) / np.pi * 180

    PWMcos_native_contact_list = []

    for i_pro in range(cg_pro_p_num):
        coor_pro_p_tmp = cg_pro_coors[i_pro]
        coors_pro_r_tmp = aagroup_pro_coors[i_pro]
        # print(i_pro, coor_pro_p_tmp, len(coors_pro_r_tmp))
        i_pro_N, i_pro_C = i_pro - 1, i_pro + 1
        if i_pro_N < 0 or cg_pro_p_segid[i_pro_N] != cg_pro_p_segid[i_pro]:
            coor_pro_p_N = coor_pro_p_tmp
        else:
            coor_pro_p_N = cg_pro_coors[i_pro_N]
        if i_pro_C >= cg_pro_p_num or cg_pro_p_segid[i_pro_C] != cg_pro_p_segid[i_pro]:
            coor_pro_p_C = coor_pro_p_tmp
        else:
            coor_pro_p_C = cg_pro_coors[i_pro_C]
        # ========= contact calculating loop =========
        for j in range(2):
            cg_dna_p_num = cg_dna_a_p_num if j == 0 else cg_dna_b_p_num
            cg_dna_coors = cg_dna_a_coors if j == 0 else cg_dna_b_coors
            aagroup_dna_coors = aagroup_dna_a_coors if j == 0 else aagroup_dna_b_coors
            cg_dna_r_ID = cg_dna_a_r_ID if j == 0 else cg_dna_b_r_ID
            cg_dna_p_ID = cg_dna_a_p_ID if j == 0 else cg_dna_b_p_ID
            for i_dna in range(1, cg_dna_p_num, 3):
                coor_dna_p_tmp = cg_dna_coors[i_dna]
                coors_dna_r_tmp = aagroup_dna_coors[i_dna]
                aa_dist_min = min_distance_residues(coors_pro_r_tmp, coors_dna_r_tmp)
                if aa_dist_min < PRO_DNA_AA_CUTOFF:
                    i_dna_5, i_dna_3 = i_dna - 3, i_dna + 3
                    if i_dna_5 < 0 or i_dna_3 > cg_dna_p_num:
                        continue
                    coor_dna_p_5 = cg_dna_coors[i_dna_5]
                    coor_dna_p_3 = cg_dna_coors[i_dna_3]
                    coor_dna_p_S = cg_dna_coors[i_dna - 1]
                    # --------- calculating r, thetas... ---------
                    v0 = coor_pro_p_tmp - coor_dna_p_tmp
                    v1 = coor_dna_p_S   - coor_dna_p_tmp
                    v2 = coor_dna_p_3   - coor_dna_p_5
                    v3 = coor_pro_p_N   - coor_pro_p_C
                    r0 = np.linalg.norm(v0)
                    theta1 = vec_angle(v0, v1)
                    theta2 = vec_angle(v0, v2)
                    theta3 = vec_angle(v0, v3)
                    PWMcos_native_contact_list.append((i_pro + 1,
                                                       cg_pro_r_ID[i_pro],
                                                       cg_dna_r_ID[i_dna],
                                                       r0, theta1, theta2, theta3,
                                                       cg_pro_p_ID[i_pro]))
    # Temporary output
    if flag_native_contact_list:
        contact_pair_output_head = "# I_AA   i_pro   i_dna         r0    theta_1    theta_2    theta_3 "
        contact_pair_output_line = "{pwmcos[0]:>6d}  {pwmcos[1]:>6d}  {pwmcos[2]:>6d} {pwmcos[3]:>10.4f} {pwmcos[4]:>10.3f} {pwmcos[5]:>10.3f} {pwmcos[6]:>10.3f} "
        print(contact_pair_output_head)
        for interaction_pair in PWMcos_native_contact_list:
            print(contact_pair_output_line.format(pwmcos=interaction_pair))


    # =======================================
    # Read in PFM and do energy decomposition
    # =======================================

    def read_pfm_to_pwm(filename):
        pfm = {}
        with open(filename, 'r') as pfm_fin:
            for line in pfm_fin:
                words = line.split()
                if len(words) < 1:
                    continue
                w0 = words[0]
                if w0 in ['A', 'C', 'G', 'T']:
                    local_list = []
                    for pfm_value in words[1:]:
                        v = float(pfm_value)
                        local_list.append(v)
                    pfm[w0] = local_list[:]
                elif w0 in {"CHAIN_A", "CHAIN_B"}:
                    local_list = []
                    for dna_res_id in words[1:]:
                        k = int(dna_res_id)
                        local_list.append(k)
                    pfm[w0] = local_list[:]
        ffp_pfm_T = np.array([pfm['A'], pfm['C'], pfm['G'], pfm['T']])
        ffp_pfm_col_sum = np.sum(ffp_pfm_T, axis=0)
        ffp_ppm_T = ffp_pfm_T / ffp_pfm_col_sum
        ffp_pwm_T_0 = -np.log(ffp_ppm_T)
        ffp_pwm_col_sum = np.sum(ffp_pwm_T_0, axis=0)
        ffp_pwm_T = ffp_pwm_T_0 - ffp_pwm_col_sum / 4
        return (ffp_pwm_T.T, pfm['CHAIN_A'], pfm['CHAIN_B'])

    pwm, dna_pwm_a_id, dna_pwm_b_id = read_pfm_to_pwm(PFM_name)

    # =========================
    # Output CafeMol ninfo file
    # =========================
    def write_ninfo():
        pwm_len = len(dna_pwm_a_id)
        ip_count = np.zeros((pwm_len, 1), dtype=np.int) 

        ninfo_name = protein_name + "_PWMcos.ninfo"
        ninfo_file = open(ninfo_name, 'w')

        ninfo_output_head = "<<<< pdpwm \n"
        ninfo_output_tail = ">>>> \n"
        ninfo_output_comm = "**     ID iunit  nil    imp        r_0    theta_1    theta_2    theta_3          ene_A          ene_C          ene_G          ene_T    gamma     eps' \n"
        ninfo_output_line = "pdpwm {pwmcosid:>3d}     3  999 {pwmcos[0]:>6d} {pwmcos[3]:>10.4f} {pwmcos[4]:>10.3f} {pwmcos[5]:>10.3f} {pwmcos[6]:>10.3f} {pwm[0]:>14.6f} {pwm[1]:>14.6f} {pwm[2]:>14.6f} {pwm[3]:>14.6f} {g:>8.4f} {e:>8.4f} \n"

        contact_pair_to_pwm = []
        for nat_contact in PWMcos_native_contact_list:
            i_dna = nat_contact[2]
            if i_dna in dna_pwm_a_id:
                contact_pair_to_pwm.append((dna_pwm_a_id.index(i_dna), 1))
                ip_count[dna_pwm_a_id.index(i_dna)] += 1
            elif i_dna in dna_pwm_b_id:
                contact_pair_to_pwm.append((dna_pwm_b_id.index(i_dna), -1))
                ip_count[dna_pwm_b_id.index(i_dna)] += 1
            else:
                print(" Index error in CHAIN_id s!   Please Check!")
                return
        for i, k in enumerate(ip_count):
            if k == 0:
                ip_count[i] = -1
        pwm_decomposed = pwm / ip_count

        ninfo_file.write(ninfo_output_head)
        ninfo_file.write(ninfo_output_comm)
        for i, nat_contact in enumerate(PWMcos_native_contact_list):
            pwm_i, pwm_v = contact_pair_to_pwm[i][0], contact_pair_to_pwm[i][1]
            pwm_line = pwm_decomposed[pwm_i][::pwm_v]
            ninfo_file.write(ninfo_output_line.format(pwmcosid=i+1, pwmcos=nat_contact, pwm=pwm_line, g=gamma, e=epsilon_prime))
        ninfo_file.write(ninfo_output_tail)
        ninfo_file.close()


    # =========================
    # Output PINANG ffp file
    # =========================
    def write_ffp():
        pwm_len = len(dna_pwm_a_id)
        ip_count = np.zeros((pwm_len, 1), dtype=np.int) 

        ffp_name = protein_name + "_PWMcos.ffp"
        ffp_file = open(ffp_name, 'w')

        ffp_output_head = "[ PWMcos ] {0:6d} \n"
        ffp_output_tail = " \n"
        ffp_output_line = "{pwmcos[7]:>6d} {pwmcos[3]:>8.4f} {pwmcos[4]:>8.3f} {pwmcos[5]:>8.3f} {pwmcos[6]:>8.3f} {pwm[0]:>9.4f} {pwm[1]:>9.4f} {pwm[2]:>9.4f} {pwm[3]:>9.4f} {g:>9.4f} {e:>9.4f}    1.0  10.0 \n"

        contact_pair_to_pwm = []
        for nat_contact in PWMcos_native_contact_list:
            i_dna = nat_contact[2]
            if i_dna in dna_pwm_a_id:
                contact_pair_to_pwm.append((dna_pwm_a_id.index(i_dna), 1))
                ip_count[dna_pwm_a_id.index(i_dna)] += 1
            elif i_dna in dna_pwm_b_id:
                contact_pair_to_pwm.append((dna_pwm_b_id.index(i_dna), -1))
                ip_count[dna_pwm_b_id.index(i_dna)] += 1
            else:
                print(" Index error in CHAIN_id s!   Please Check!")
                return
        for i, k in enumerate(ip_count):
            if k == 0:
                ip_count[i] = -1
        pwm_decomposed = pwm / ip_count

        ffp_file.write(ffp_output_head.format(len(PWMcos_native_contact_list)))
        for i, nat_contact in enumerate(PWMcos_native_contact_list):
            pwm_i, pwm_v = contact_pair_to_pwm[i][0], contact_pair_to_pwm[i][1]
            pwm_line = pwm_decomposed[pwm_i][::pwm_v]
            ffp_file.write(ffp_output_line.format(pwmcos=nat_contact, pwm=pwm_line, g=gamma, e=epsilon_prime))
        ffp_file.write(ffp_output_tail)
        ffp_file.close()

    write_ffp()
    write_ninfo()


if __name__ == '__main__':
    import argparse
    def parse_arguments():
        parser = argparse.ArgumentParser(description='Generate CafeMol ninfo file from PDB and PFM.')
        parser.add_argument('pdb', type=str, help="PDB file name.")
        parser.add_argument('pfm', type=str, help="PFM (position frequency matrix) file name.")
        parser.add_argument('-t', '--psf', help="Output psf file.", action="store_true")
        parser.add_argument('-n', '--contact', help="Print native contact list.", action="store_true")
        parser.add_argument('-c', '--cutoff', type=float, default=4.0, help="Specify cutoff in native contact definition.")
        parser.add_argument('-g', '--gamma', type=float, default=1.0, help="Specify Gamma, the energy scaling factor.")
        parser.add_argument('-e', '--epsilon', type=float, default=0.0, help="Specify Epsilon', the energy shift factor.")
        return parser.parse_args()
    args = parse_arguments()
    main(args.pdb, args.pfm, args.psf, args.contact, args.gamma, args.epsilon, args.cutoff)

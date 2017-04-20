#!/usr/bin/env python3

from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def read_sequence(fasta_file_name):
    """Read DNA sequence from .fasta file.
    Keyword Arguments:
    fasta_file -- sequence file in the format of .fasta
    """
    fasta_sequences = SeqIO.parse(open(fasta_file_name),'fasta')
    dna_fasta = list(fasta_sequences)[0]
    return dna_fasta

def read_pwm(pwm_file_name):
    """Read position weight matrix from pwm_file.
    Keyword Arguments:
    pwm_file_name -- 
    """
    pwm = {}
    with open(pwm_file_name, 'r') as pwm_fin:
        for line in pwm_fin:
            words = line.split()
            base_type = line[0]
            enescore = []
            for v in words[1:]:
                enescore.append(float(v))
            pwm[base_type] = enescore[:]
    return pwm.copy()

def main(fasta_file_name, pwm_file_name):
    # ------------------------------ read fasta ------------------------------
    dna_fasta = read_sequence(fasta_file_name)
    dna_name, dna_sequence = dna_fasta.id, dna_fasta.seq
    reverse_dna_sequence = dna_sequence.reverse_complement()
    # print(dna_name, dna_sequence)
    # print(dna_sequence.reverse_complement())
    dna_len = len(dna_sequence)

    # ------------------------------ read pwm ------------------------------
    pwm = read_pwm(pwm_file_name)
    pwm_len = len(pwm['A'])

    # ------------------------------ pwm score calculation ------------------------------
    def calculate_pwm_energy(seq_piece, pwm0):
        score = 0
        for i, b in enumerate(seq_piece):
            s = pwm0[b][i]
            score += s * 0.593
        return score
    
    # ------------------------------ simple test ------------------------------
    shifting_score = []
    shifting_score_reverse = []
    if dna_len < pwm_len:
        print("ERROR: DNA length < PWM length!!! STOP here.")
        return 1
    piece_start, piece_end = 0, pwm_len
    while piece_end <= dna_len:
        dna_piece = dna_sequence[piece_start : piece_end]
        reverse_dna_piece = dna_piece.reverse_complement()
        piece_score = calculate_pwm_energy(dna_piece, pwm)
        reverse_piece_score = calculate_pwm_energy(reverse_dna_piece, pwm)
        # print(piece_start, dna_piece, reverse_dna_piece, piece_score, reverse_piece_score)
        shifting_score.append(piece_score)
        shifting_score_reverse.append(reverse_piece_score)
        piece_start += 1
        piece_end += 1

    # ------------------------------ plotting ------------------------------
    fig, axes = plt.subplots(1, 1, figsize=(12, 8))
    X = [i for i in range(dna_len - pwm_len + 1)]
    axes.plot(X, shifting_score, 'r')
    axes.plot(X, shifting_score_reverse, 'g')
    axes.set_xlim(0, dna_len - pwm_len)
    # axes.set_ylim(-12, 10)
    # axes.set_xticks()
    # axes.set_yticks()
    # axes.text(x=0, y=0, s=str(dna_sequence), family='monospace', fontsize=4.6)
    # axes.text(x=0, y=0.5, s=str(reverse_dna_sequence), family='monospace', fontsize=4.6)

    plt.show()
    

if __name__ == '__main__':
    import sys
    def print_usage():
        usage = sys.argv[0] + " seq.fasta tf.pwm"
    fasta_file_name = sys.argv[1]
    pwm_file_name = sys.argv[2]
    main(fasta_file_name, pwm_file_name)

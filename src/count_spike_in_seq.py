#!/usr/bin/env python2.7
import fastqutil
import sys
import collections

# The length to skip before checking spike in sequences. 
PFX_SKIP_LEN = 15

# Spike-in sequence list
#SS_LIST = ['GAACAGTATTTGGTATCTGCGCTCTGCACG',
#           'AACAGTATTTGGTATCTGCGCTCTGCACGA',
#           'ACAGTATTTGGTATCTGCGCTCTGCACGAT', 
#           'CAGTATTTGGTATCTGCGCTCTGCACGATG',
#           'AGTATTTGGTATCTGCGCTCTGCACGATGG']

# Extended spike-in sequence list
SS_LIST = ['GAAGAACAGTATTTGGTATCTGCGCTCTGC',
           'AAGAACAGTATTTGGTATCTGCGCTCTGCA',
           'AGAACAGTATTTGGTATCTGCGCTCTGCAC',
           'GAACAGTATTTGGTATCTGCGCTCTGCACG',
           'AACAGTATTTGGTATCTGCGCTCTGCACGA',
           'ACAGTATTTGGTATCTGCGCTCTGCACGAT',
           'CAGTATTTGGTATCTGCGCTCTGCACGATG',
           'AGTATTTGGTATCTGCGCTCTGCACGATGG',
           'GTATTTGGTATCTGCGCTCTGCACGATGGG',
           'TATTTGGTATCTGCGCTCTGCACGATGGGT',
           'ATTTGGTATCTGCGCTCTGCACGATGGGTT',
           'TTTGGTATCTGCGCTCTGCACGATGGGTTA']

# Get spike in sequence length from spike in sequence list.
# Assert that all spike in sequences have the same length.
def get_ss_len(ss_list):
    assert len(ss_list) >= 1
    ss_len = len(ss_list[0])

    for ss in ss_list:
        assert len(ss) == ss_len, "All spike-in sequences should have the same length."

    return ss_len

# Get sequence substring for spike-in sequence checking.
def get_ss_check_seq(seq, pfx_skip_len, ss_len):
    return seq[pfx_skip_len:pfx_skip_len + ss_len]

# Use a dictionary to count the number of spike in sequences
def count_spike_in_seq(ifn, ofn):
    ss_len = get_ss_len(SS_LIST)
    
    ss_cnt_dict = {}
    for ss in SS_LIST:
        ss_cnt_dict[ss] = 0

    for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(ifn):
        ss_check_seq = get_ss_check_seq(seq, PFX_SKIP_LEN, ss_len)
        if ss_check_seq in ss_cnt_dict:
            ss_cnt_dict[ss_check_seq] += 1

    with open(ofn, 'w') as ofile:
        ofile.write('Spike-in Sequence\tCount\n')
        for ss in SS_LIST:
            ofile.write(ss + '\t' + str(ss_cnt_dict[ss]) + '\n')

    return
    
def main():
    argv = sys.argv
    if len(argv) != 3:
        sys.stderr.write("Usage:\n%s <input fastq file> <output stats file>\n" % argv[0])
        sys.exit(-1)

    ifn = argv[1]
    ofn = argv[2]
    count_spike_in_seq(ifn, ofn)

if __name__ == '__main__':
    main()
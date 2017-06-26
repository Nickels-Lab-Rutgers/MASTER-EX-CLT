#!/usr/bin/env python2.7

import sys
import itertools
import rnautil

BEG_SEQ = "GCGGCCGCAGGCTTGACACTTTATGCTTCGGCTCGTATAATGTG"

def ex_count_parsed_rna(dna_rd_len, dna_rd_start_ind, i_fn, mm_fn, o_fn):
    # DNA region includes extra sequence downstream
    # dna_prmt_seq_len contains the suffix length
    dna_prmt_seq_len = rnautil.get_dna_prmt_seq_len(i_fn)
    max_rna_prmt_seq_len = dna_prmt_seq_len + 5
    count_dict = {}
    # {N10 : {'TNM' : {len : count, ...}, 
    #         'TNA' : {len : count, ...}, 
    #         ...},
    # ...}
    for prmt_seq_tup in itertools.product("ACGT", repeat = dna_rd_len):
        cnt_type_list = ['TNM', 'TNA', 'M', 'A']
        count_dict[''.join(prmt_seq_tup)] = dict(zip(cnt_type_list, 
            [dict(zip(range(1, max_rna_prmt_seq_len + 1), 
                      [0] * max_rna_prmt_seq_len)) for i in xrange(len(cnt_type_list))]))

    i_file = open(i_fn)

    for line in i_file:
        fields = line.strip().split()
        assert len(fields) == 6

        rna = fields[0]
        dna = fields[1]
        dna_rd_seq = dna[dna_rd_start_ind : dna_rd_start_ind + dna_rd_len]
        raw_cnt = int(fields[3])
        tn_cnt = int(fields[4])
        rna_len = len(rna)
        assert rna_len > 0

        # Only count RNA with promoter seq len [1, DNA prmt seq len + 5]
        if rna_len <= max_rna_prmt_seq_len:
            count_dict[dna_rd_seq]['TNA'][rna_len] += tn_cnt
            count_dict[dna_rd_seq]['A'][rna_len] += raw_cnt
            # rd_match = rnautil.check_match(dna, rna)
            dna_rna_match = rnautil.check_match(BEG_SEQ + dna, rna)
            if dna_rna_match:
                count_dict[dna_rd_seq]['TNM'][rna_len] += tn_cnt
                count_dict[dna_rd_seq]['M'][rna_len] += raw_cnt

    i_file.close()

    o_file = open(o_fn, 'w')
    mm_file = open(mm_fn, 'w')

    for prmt_seq in sorted(count_dict.keys()):
        o_file.write(prmt_seq)
        for cnt_type in ['TNM', 'TNA', 'M', 'A']:
            for rna_len in sorted(count_dict[prmt_seq][cnt_type].keys(), reverse = True):
                o_file.write('\t' + str(count_dict[prmt_seq][cnt_type][rna_len]))
        o_file.write('\n')

        # Calculate mismatch percent. sum(TNM[1:dna_prmt_len]) / sum(TNA[1:max_rna_len])
        # Tx before random region is counted as mismatch
        tnm_prmt_sum = 0
        tna_prmt_sum = 0
        for rna_len in xrange(1, dna_prmt_seq_len + 1):
            tnm_prmt_sum += count_dict[prmt_seq]['TNM'][rna_len]

        tna_prmt_sum = sum(count_dict[prmt_seq]['TNA'].values())

        assert tnm_prmt_sum <= tna_prmt_sum

        if tna_prmt_sum > 0:
            mm_file.write("%s\t%f\n" % (prmt_seq, 
                (tna_prmt_sum - tnm_prmt_sum) / float(tna_prmt_sum) * 100))
        else:
            mm_file.write(prmt_seq + '\n')

    o_file.close()
    mm_file.close()
    return

def main():
    argv = sys.argv
    if len(argv) != 6:
        sys.stderr.write("Usage:\n%s [DNA random region length] [DNA random region start index (0-based)] \
<input RNA parsed file> <output mismatch percentage file> \
<output RNA count file> \nSchema for output is Promoter Seq, TNM, TNA, M, and A. \n\
Length of counted reads ranges from 1 to (DNA promoter region length + 5)\n" % argv[0])
        return -1

    dna_rd_len = int(argv[1])
    dna_rd_start_ind = int(argv[2])
    if dna_rd_start_ind < 0:
        sys.stderr.write("[DNA random region start index (0-based)] should >= 0")
        return -2

    i_fn = argv[3]
    mm_fn = argv[4]
    o_fn = argv[5]
    ex_count_parsed_rna(dna_rd_len, dna_rd_start_ind, i_fn, mm_fn, o_fn)


if __name__ == "__main__":
    main()

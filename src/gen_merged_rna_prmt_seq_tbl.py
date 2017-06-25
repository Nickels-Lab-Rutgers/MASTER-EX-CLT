#!/usr/bin/env python2.7

import sys
import rnautil


def get_mp_dict(mmfn):
    mp_dict = {}
    # mp_file format: DNA prmt seq \t mismatch percentage.
    # If TNA sum is 0, there is not \t and mismatch percentage. 
    with open(mmfn) as mp_file:
        for line in mp_file:
            fields = line.strip().split()
            if len(fields) > 1:
                mp_dict[fields[0]] = fields[1]
            else:
                mp_dict[fields[0]] = 'NonApplicable'

    return mp_dict
    
def gen_rna_prmt_summary(mmfn, i_fn, o_fnm, o_fnl, o_fna):
    mp_dict = get_mp_dict(mmfn)

    dna_prmt_seq_len = rnautil.get_dna_prmt_seq_len(i_fn)
    max_rna_prmt_seq_len = dna_prmt_seq_len + 5

    # { (dna_prmt, rna_prmt) : {'raw_cnt' : 0, 'tn_cnt' : 0}, ...}
    rna_prmt_dict = {}

    i_file = open(i_fn)
    o_filem = open(o_fnm, 'w')
    o_filel = open(o_fnl, 'w')
    o_filea = open(o_fna, 'w')

    # parsed_start_pos and parsed_match_bit are obsolete fields.
    for (rna, dna, parsed_start_pos, raw_cnt, 
         tn_cnt, parsed_match_bit) in rnautil.iterate_rna_parsed_file(i_fn):
        assert len(dna) == dna_prmt_seq_len

        if len(rna) > max_rna_prmt_seq_len:
            continue

        if (dna, rna) not in rna_prmt_dict:
            rna_prmt_dict[(dna, rna)] = {'raw_cnt' : int(raw_cnt), 
                                         'tn_cnt' : int(tn_cnt)}
        else:
            rna_prmt_dict[(dna, rna)]['raw_cnt'] += int(raw_cnt)
            rna_prmt_dict[(dna, rna)]['tn_cnt'] += int(tn_cnt)

    for (dna, rna) in sorted(rna_prmt_dict.keys()):
        tn_cnt = rna_prmt_dict[(dna, rna)]['tn_cnt']

        tx_start_pos = dna_prmt_seq_len - len(rna) + 1
        seq_match = rnautil.check_match(dna, rna)

        if seq_match:
            mmark = 'match'
        else:
            mmark = 'mismatch'

        if tx_start_pos <= 0:
            lmark = 'L'
        else:
            lmark = 'N' + str(dna_prmt_seq_len)

        mmpct = mp_dict[dna]

        rna_prmt_line = '\t'.join((dna, rna, mmark, lmark, str(tx_start_pos), 
                                   str(tn_cnt), mmpct)) + '\n'
        
        o_filea.write(rna_prmt_line)
        if not seq_match:
            o_filem.write(rna_prmt_line)
        if tx_start_pos <= 0:
            o_filel.write(rna_prmt_line)


    o_filem.close()
    o_filel.close()
    o_filea.close()
    return


def main():
    argv = sys.argv
    if len(argv) != 6:
        print "Usage:\n%s <input RNA mismatch percentages file> <input RNA parsed file> \
<output mismatch RNA promoter region summary file> <output L RNA promoter region summary file> \
<output all RNA promoter region summary file>\n\
Schema for output is DNA Tpl, RNA Seq, Start Pos, Match?, TN Count, Count." %argv[0]
        sys.exit(2)
    mmfn = argv[1]
    i_fn = argv[2]
    o_fnm = argv[3]
    o_fnl = argv[4]
    o_fna = argv[5]
    gen_rna_prmt_summary(mmfn, i_fn, o_fnm, o_fnl, o_fna)


if __name__ == "__main__":
    main()
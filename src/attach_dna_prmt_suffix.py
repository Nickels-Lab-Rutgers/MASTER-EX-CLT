#!/usr/bin/env python2.7
import sys
import itertools

# Parsed DNA format:
# barcode promoter_region count
def is_valid_dna_seq(seq):
    for nt in seq:
        if nt not in 'ACGT':
            return False
    return True

def attach_dna_prmt_suffix(ex_prmt_suffix_seq, i_fn, o_fn):
    i_file = open(i_fn, 'r')
    o_file = open(o_fn, 'w')
    
    for line in i_file:
        fields = line.strip().split()
        brcd_seq = fields[0]
        prmt_seq = fields[1]
        cnt = fields[2]

        ofile.write('\t'.join((brcd_seq, prmt_seq + ex_prmt_suffix_seq, cnt)) + '\n')

    i_file.close()
    o_file.close()

    return

def main():
    argv = sys.argv
    if len(argv) != 4:
        sys.stderr.write("Usage:\n%s [suffix sequence] <input DNA template file> \
<output DNA template file with extra suffix>\n" % argv[0])
        return(-1)

    ex_prmt_suffix_seq = argv[1]
    if not is_valid_dna_seq(ex_prmt_suffix_seq):
        sys.stderr.write("Invalid suffix sequence: %s\nSuffix sequence can only contain ACGT. \n" % ex_prmt_suffix_seq)
        return(-2)
    
    i_fn = argv[2]
    o_fn = argv[3]
    attach_dna_prmt_suffix(ex_prmt_suffix_seq, i_fn, o_fn)


if __name__ == '__main__':
    main()
#!/usr/bin/env python2.7
import sys
import itertools

# Parsed DNA format:
# barcode promoter_region count
def get_promoter_region_len(i_fn):
    i_file = open(i_fn)
    line = i_file.readline()
    i_file.close()

    fields = line.strip().split()
    assert len(fields) == 3

    promoter_region = fields[1]
    return(len(promoter_region))
    

def gen_dna_count(i_fn, o_fn):
    prmt_region_len = get_promoter_region_len(i_fn)
    i_file = open(i_fn)
    cnt_dict = {}
    
    for prmt_seq_tup in itertools.product("ACGT", repeat = prmt_region_len):
        cnt_dict[''.join(prmt_seq_tup)] = 0

    for line in i_file:
        fields = line.strip().split()
        prmt_seq = fields[1]
        cnt = int(fields[2])
        # Implicitly, all promoter sequences must have the same length. 
        cnt_dict[prmt_seq] += cnt

    i_file.close()

    o_file = open(o_fn, 'w')
    for tag in sorted(cnt_dict.keys()):
        o_file.write(tag + '\t' + str(cnt_dict[tag]) + '\n')

    o_file.close()

    return

def main():
    argv = sys.argv
    if len(argv) != 3:
        print "Usage:\n%s <input DNA template file> <output DNA count file>" % argv[0]
        return(-1)
    i_fn = argv[1]
    o_fn = argv[2]
    gen_dna_count(i_fn, o_fn)


if __name__ == '__main__':
    main()
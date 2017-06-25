#!/usr/bin/env python2.7
import fastqutil
import sys
import collections

def count_seq_pfx(num_pfx_bp, ifn, ofn):
    num_total_reads = 0
    seq_pfx_list = []
    for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(ifn):
        num_total_reads += 1
        if len(seq) < num_pfx_bp:
            raise ValueError('Read %s has less base pairs than the number of \
base pairs for counting.' % seqid)

        seq_pfx_list.append(seq[:num_pfx_bp])

    seq_pfx_cnter = collections.Counter(seq_pfx_list)

    pfx_cnt_list = [(pfx, seq_pfx_cnter[pfx]) for pfx in seq_pfx_cnter]
    sorted_pfx_cnt_list = sorted(pfx_cnt_list, key = lambda tup: -tup[1])

    cnt_sum = sum([seq_pfx_cnter[pfx] for pfx in seq_pfx_cnter])
    assert cnt_sum == num_total_reads

    with open(ofn, 'w') as ofile:
        ofile.write('Sequence Prefix\tCount\tPercentage\n')
        for pfx_cnt_tup in sorted_pfx_cnt_list:
            ofile.write("%s\t%d\t%f\n" % (pfx_cnt_tup[0], pfx_cnt_tup[1], pfx_cnt_tup[1] / float(num_total_reads) * 100))



def main():
    argv = sys.argv
    if len(argv) != 4:
        sys.stderr.write("Usage:\n%s [Number of prefix bp to count] \
<input fastq file> <output stats file>\n" % argv[0])
        sys.exit(-1)

    num_pfx_bp = int(argv[1])
    ifn = argv[2]
    ofn = argv[3]
    count_seq_pfx(num_pfx_bp, ifn, ofn)



if __name__ == '__main__':
    main()

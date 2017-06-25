#!/usr/bin/env python2.7
import sys
import fastqutil

def is_valid_pfx_sep(pfx_sep):
    for c in pfx_sep:
        if c not in 'ACGT':
            return False
    return True

def prefix_grep_fastq(pfx_sep, ifn, ofn):
    assert is_valid_pfx_sep(pfx_sep), 'Prefix separator must only contain ACGT'

    o_fastq_file = open(ofn, 'w')
    for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(ifn):
        if seq[:len(pfx_sep)] == pfx_sep:
            o_fastq_file.write('\n'.join((seqid, seq, rseqid, qscore)) + '\n')
    o_fastq_file.close()
    return

def main():
    argv = sys.argv
    if len(argv) != 4:
        sys.stderr.write('Usage:\n%s <prefix separator> <input fastq file> \
<output fastq file>\n' % argv[0])
        return -2

    pfx_sep = argv[1]
    ifn = argv[2]
    ofn = argv[3]

    prefix_grep_fastq(pfx_sep, ifn, ofn)

if __name__ == "__main__":
    main()

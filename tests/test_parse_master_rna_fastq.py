#!/usr/bin/env python2.7
import env
import unittest
from src import parse_master_rna_fastq

def construct_rna_fastq_rec(dtag, brcd, prmt, count):
    seq = (dtag + prmt + parse_master_rna_fastq.MID_SEQ + brcd +
           parse_master_rna_fastq.END_SEQ)
    qscore = chr(33 + 33) * len(seq)
    rec = "@seqid\n%s\n+rseqid\n%s\n" % (seq, qscore)
    return (rec * count)


class TestFastqFileParser(unittest.TestCase):
    def test_linenum(self):
        self.assertEqual(
            parse_master_rna_fastq.parse_rna_fastq_files('!', 15, 7, 15, "data/dnatpl.txt",
                "data/testout.txt", "data/teststat.txt", ["data/linenum.txt"]),
            3)

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_fastq_files('!', 15, 7, 15, "data/dnatpl.txt",
                "data/testout.txt", "data/teststat.txt", ["data/empty_line.txt"]), 
            2)

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_fastq_files('!', 15, 7, 15, "data/dnatpl.txt",
                "data/testout.txt", "data/teststat.txt", ["data/linenum.txt", "data/empty_line.txt"]), 
            5)

    def test_bad_fastq_file(self):
        with self.assertRaises(ValueError):
            parse_master_rna_fastq.parse_rna_fastq_files('!', 15, 7, 15, "data/dnatpl.txt",
                "data/testout.txt", "data/teststat.txt", ["data/line_num_not_multiple_of_4.txt"])

    def test_result(self):
        dna_tpl_list = [('CGTACGTACGTACGT', 'ACGTACG', 10),
                        ('AAAAAAAAAAAAAAA', 'CCCCCCC', 10)]
        with open("data/dnatpl.txt", 'w') as ofile:
            ofile.write('\n'.join(['\t'.join(map(str, tup)) for tup in dna_tpl_list]) + '\n')

        rna_seq_list = sorted([('TTTTTTTTTTTTTTT', 'CC', 'AAAAAAAAAAAAAAA', 10),
                               ('TTTTTTTTTTTTTTC', 'CC', 'AAAAAAAAAAAAAAA', 10),
                               ('TTTTTTTTTTTTTTG', 'CC', 'AAAAAAAAAAAAAAA', 10),
                               ('TTTTTTTTTTTTTTG', 'GG', 'AAAAAAAAAAAAAAA', 30),
                               ('TTTTTTTTTTTTTTG', 'ACGTACG', 'CGTACGTACGTACGT', 20),
                               ('TTTTTTTTTTTTTTA', 'CACGTACG', 'CGTACGTACGTACGT', 1),
                               ('TTTTTTTTTTTTTTA', 'CAGGTACG', 'CGTACGTACGTACGT', 2),
                               ('TTTTTTTTTTTTTAA', 'CCACGTACG', 'CGTACGTACGTACGT', 3),
                               ('TTTTTTTTTTTTTAA', 'CCCAAGTACG', 'CGTACGTACGTACGT', 5)
                              ],
                              key = lambda tup: tup[2])
        rna_result_list = [('CC', 'CCCCCCC', 6, 30, 3, 1), 
                           ('GG', 'CCCCCCC', 6, 30, 1, 0),
                           ('ACGTACG', 'ACGTACG', 1, 20, 1, 1),
                           ('CACGTACG', 'ACGTACG', 0, 1, 1, 1),
                           ('CAGGTACG', 'ACGTACG', 0, 2, 1, 0),
                           ('CCACGTACG', 'ACGTACG', -1, 3, 1, 1),
                           ('CCCAAGTACG', 'ACGTACG', -2, 5, 1, 0),
                           ]

        ffile = open('data/rna1.fastq', 'w')
        for dtagseq, prmtseq, brcdseq, count in rna_seq_list:
            ffile.write(construct_rna_fastq_rec(dtagseq, brcdseq, prmtseq, count))
        ffile.close()

        parse_master_rna_fastq.parse_rna_fastq_files('!', 15, 7, 15, "data/dnatpl.txt",
            "data/rnap1test.txt", "data/rnas1test.txt", ["data/rna1.fastq"])
        with open('data/rnap1test.txt', 'r') as ifile:
            r1testlist = sorted([line.strip().split() for line in ifile])

        r1reflist = sorted([map(str, tup) for tup in rna_result_list])

        self.assertEqual(r1testlist, r1reflist)


class TestRNASeqParser(unittest.TestCase):
    tpl_seq = ('C' * 15 + 'AAACGTACG' + parse_master_rna_fastq.MID_SEQ + 
        'CGTACGTACGTACGT' + parse_master_rna_fastq.END_SEQ)
    padded_tpl_seq = tpl_seq + 'AGGTGGAATTCTCGGGTGCCAAGG'
    def test_correct_seq(self):
        self.assertEqual(parse_master_rna_fastq.parse_rna_seq(self.tpl_seq, 
                chr(50) * len(self.tpl_seq), chr(33), 15, 7, 15),
            ('C' * 15, 'AAACGTACG', 'CGTACGTACGTACGT'))
        self.assertEqual(parse_master_rna_fastq.parse_rna_seq(self.padded_tpl_seq, 
                chr(50) * len(self.padded_tpl_seq), chr(33), 15, 7, 15),
            ('C' * 15, 'AAACGTACG', 'CGTACGTACGTACGT'))

    def test_qual_failure(self):
        self.assertEqual(parse_master_rna_fastq.parse_rna_seq(self.tpl_seq, 
                chr(39) + chr(50) * (len(self.tpl_seq) - 1), chr(39), 15, 7, 15),
            ('C' * 15, 'AAACGTACG', 'CGTACGTACGTACGT'))
        self.assertEqual(parse_master_rna_fastq.parse_rna_seq(self.tpl_seq, 
                chr(38) + chr(50) * (len(self.tpl_seq) - 1), chr(39), 15, 7, 15),
            'QUAL_FAILED')
        self.assertEqual(parse_master_rna_fastq.parse_rna_seq(
            self.tpl_seq.replace('ACGTACG', 'ACNTACG'), 
                chr(39) + chr(50) * (len(self.tpl_seq) - 1), chr(39), 15, 7, 15),
            'QUAL_FAILED')


    def test_structure_failure(self):
        print '\ntpl seq: ' + self.tpl_seq
        seq_mid_mm = ('C' * 15 + 'AAACGTACG' + 'C' * 5 + 
            parse_master_rna_fastq.MID_SEQ[5:] + 'CGTACGTACGTACGT' + 
            parse_master_rna_fastq.END_SEQ)
        print 'mid mm : ' + seq_mid_mm

        seq_end_mm = ('C' * 15 + 'AAACGTACG' + 
            parse_master_rna_fastq.MID_SEQ + 'CGTACGTACGTACGT' + 'A' + 
            parse_master_rna_fastq.END_SEQ[1:])
        print 'end mm : ' + seq_end_mm

        seq_short = seq_end_mm[:-1]
        print 'short  : ' + seq_short

        seq_p_short = ('C' * 15 + 'ACGTAC' + 
            parse_master_rna_fastq.MID_SEQ + 'CGTACGTACGTACGT' +
            parse_master_rna_fastq.END_SEQ + 'AAA')
        print 'pshort : ' + seq_p_short

        seq_b_short = ('C' * 15 + 'ACGTACG' + 
            parse_master_rna_fastq.MID_SEQ + 'CGTACGTACGTACG' + 'A' + 
            parse_master_rna_fastq.END_SEQ[1:] + 'AAA')
        print 'bshort : ' + seq_b_short

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_mid_mm, chr(50) * len(seq_mid_mm),
                chr(33), 15, 7, 15),
            'MID_SEQ_NOT_FOUND')

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_end_mm, chr(50) * len(seq_end_mm),
                chr(33), 15, 7, 15),
            'MID_SEQ_NOT_FOUND')

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_short, chr(50) * len(seq_short),
                chr(33), 15, 7, 15),
            'SHORT_SEQ')

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_p_short, chr(50) * len(seq_p_short),
                chr(33), 15, 7, 15),
            ('C' * 15, 'ACGTAC', 'CGTACGTACGTACGT'))

        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_b_short, chr(50) * len(seq_b_short),
                chr(33), 15, 7, 15),
            'MID_SEQ_NOT_FOUND')

        seq_d_short = ('C' * 15 + parse_master_rna_fastq.MID_SEQ + 'CGTACGTACGTACGT' +
            parse_master_rna_fastq.END_SEQ + 'AAA')
        print 'dshort : ' + seq_d_short
        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_d_short, chr(50) * len(seq_d_short),
                chr(33), 15, 7, 15),
            'MID_SEQ_NOT_FOUND')

        seq_p1 = ('C' * 15 + 'L' + parse_master_rna_fastq.MID_SEQ + 'CGTACGTACGTACGT' +
            parse_master_rna_fastq.END_SEQ + 'AAA')
        print 'seq_p1 : ' + seq_p1
        self.assertEqual(
            parse_master_rna_fastq.parse_rna_seq(seq_p1, chr(50) * len(seq_p1),
                chr(33), 15, 7, 15),
            ('C' * 15, 'L', 'CGTACGTACGTACGT'))


if __name__ == '__main__':
    unittest.main()


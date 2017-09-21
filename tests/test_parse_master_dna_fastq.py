#!/usr/bin/env python2.7
import env
import unittest
from src import parse_master_dna_fastq
# Assume that the read count cutoff in dna parser is 10. 


def construct_dna_fastq_rec(brcd, prmt, count):
    seq = (parse_master_dna_fastq.BEG_SEQ + prmt + 
           parse_master_dna_fastq.MID_SEQ + brcd + 
           parse_master_dna_fastq.END_SEQ)
    qscore = chr(33 + 33) * len(seq)
    rec = "@seqid\n%s\n+rseqid\n%s\n" % (seq, qscore)
    return (rec * count)

class TestFastqFileParser(unittest.TestCase):
    def test_linenum(self):
        self.assertEqual(
            parse_master_dna_fastq.parse_dna_fastq_files('!', 7, 15, "data/testout.txt", 
                "data/dna_teststat.txt", ["data/linenum.txt"]), 
            3)

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_fastq_files('!', 7, 15, "data/testout.txt", 
                "data/dna_teststat.txt", ["data/empty_line.txt"]), 
            2)

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_fastq_files('!', 7, 15, "data/testout.txt", 
                "data/dna_teststat.txt", ["data/linenum.txt", "data/empty_line.txt"]), 
            5)

    def test_bad_fastq_file(self):
        with self.assertRaises(ValueError):
            parse_master_dna_fastq.parse_dna_fastq_files('!', 7, 15, "data/testout.txt", 
                "data/dna_teststat.txt", ["data/line_num_not_multiple_of_4.txt"])

    def test_result(self):
        brcd_prmt_list = sorted([('CGTACGTACGTACGT', 'ACGTACG', 10),
                                 ('AAAAAAAAAAAAAAA', 'CCCCCCC', 10),
                                 ('AAAAAAAAAAAAAAC', 'CCCCCCC', 9),
                                 ('AAAAAAAAAAAAAAG', 'CCCCCCC', 9),
                                 ('AAAAAAAAAAAAAAG', 'CCCCCCG', 1),
                                 ('AAAAAAAAAAAAAAG', 'CCCCCCT', 1),
                                 ('AAAAAAAAAAAAAAT', 'ACCCCCC', 10),
                                 ('AAAAAAAAAAAAAAT', 'CCCCCCC', 3)],
                                 key = lambda tup: tup[0])
        result_list    = sorted([('CGTACGTACGTACGT', 'ACGTACG', 10),
                                 ('AAAAAAAAAAAAAAA', 'CCCCCCC', 10)],
                                 key = lambda tup: tup[0])
        ffile = open('data/dna1.fastq', 'w')
        for brcdseq, prmtseq, count in brcd_prmt_list:
            ffile.write(construct_dna_fastq_rec(brcdseq, prmtseq, count))
        ffile.close()
            
        parse_master_dna_fastq.parse_dna_fastq_files('!', 7, 15, "data/dnap1test.txt", 
            "data/dnas1test.txt", ["data/dna1.fastq"])
        with open("data/dnap1test.txt") as ifile:
            pteststr = ifile.read()

        prefstr = '\n'.join(['\t'.join(map(str, tup)) for tup in result_list]) + '\n'

        self.assertEqual(prefstr, pteststr)


class TestDNASeqParser(unittest.TestCase):
    tpl_seq = (parse_master_dna_fastq.BEG_SEQ + 'ACGTACG' + parse_master_dna_fastq.MID_SEQ + 
        'CGTACGTACGTACGT' + parse_master_dna_fastq.END_SEQ)
    padded_tpl_seq = 'GCGGCCGC' + tpl_seq + 'AGGTGGAATTCTCGGGTGCCAAGG'
    def test_correct_seq(self):
        self.assertEqual(parse_master_dna_fastq.parse_dna_seq(self.tpl_seq, 
                chr(50) * len(self.tpl_seq), chr(33), 7, 15),
            ('ACGTACG', 'CGTACGTACGTACGT'))
        self.assertEqual(parse_master_dna_fastq.parse_dna_seq(self.padded_tpl_seq, 
                chr(50) * len(self.padded_tpl_seq), chr(33), 7, 15),
            ('ACGTACG', 'CGTACGTACGTACGT'))

    def test_qual_failure(self):
        self.assertEqual(parse_master_dna_fastq.parse_dna_seq(self.tpl_seq, 
                chr(39) + chr(50) * (len(self.tpl_seq) - 1), chr(39), 7, 15),
            ('ACGTACG', 'CGTACGTACGTACGT'))
        self.assertEqual(parse_master_dna_fastq.parse_dna_seq(self.tpl_seq, 
                chr(38) + chr(50) * (len(self.tpl_seq) - 1), chr(39), 7, 15),
            'QUAL_FAILED')
        self.assertEqual(parse_master_dna_fastq.parse_dna_seq(
            self.tpl_seq.replace('ACGTACG', 'ACNTACG'), 
                chr(39) + chr(50) * (len(self.tpl_seq) - 1), chr(39), 7, 15),
            'QUAL_FAILED')


    def test_structure_failure(self):
        print '\ntpl seq: ' + self.tpl_seq
        seq_beg_mm = 'A' * 5 + self.tpl_seq[5:]
        print 'beg mm : ' + seq_beg_mm
        seq_mid_mm = (parse_master_dna_fastq.BEG_SEQ + 'ACGTACG' + 'C' * 5 + 
            parse_master_dna_fastq.MID_SEQ[5:] + 'CGTACGTACGTACGT' + 
            parse_master_dna_fastq.END_SEQ)
        print 'mid mm : ' + seq_mid_mm
        seq_end_mm = (parse_master_dna_fastq.BEG_SEQ + 'ACGTACG' + 
            parse_master_dna_fastq.MID_SEQ + 'CGTACGTACGTACGT' + 'A' + 
            parse_master_dna_fastq.END_SEQ[1:])
        print 'end mm : ' + seq_end_mm

        seq_short = seq_mid_mm[:-1]
        print 'short  : ' + seq_short
        seq_p_short = (parse_master_dna_fastq.BEG_SEQ + 'ACGTAC' + 
            parse_master_dna_fastq.MID_SEQ + 'CGTACGTACGTACGT' + 'A' + 
            parse_master_dna_fastq.END_SEQ[1:] + 'AAA')
        print 'pshort : ' + seq_p_short
        seq_b_short = (parse_master_dna_fastq.BEG_SEQ + 'ACGTACG' + 
            parse_master_dna_fastq.MID_SEQ + 'CGTACGTACGTACG' + 'A' + 
            parse_master_dna_fastq.END_SEQ[1:] + 'AAA')
        print 'bshort : ' + seq_b_short
        
        self.assertEqual(
            parse_master_dna_fastq.parse_dna_seq(seq_beg_mm, chr(50) * len(seq_beg_mm),
                chr(33), 7, 15),
            'BEG_SEQ_NOT_FOUND')

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_seq(seq_mid_mm, chr(50) * len(seq_mid_mm),
                chr(33), 7, 15),
            'MID_SEQ_NOT_FOUND')

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_seq(seq_end_mm, chr(50) * len(seq_end_mm),
                chr(33), 7, 15),
            'END_SEQ_NOT_FOUND')

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_seq(seq_short, chr(50) * len(seq_short),
                chr(33), 7, 15),
            'SHORT_SEQ')

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_seq(seq_p_short, chr(50) * len(seq_p_short),
                chr(33), 7, 15),
            'MID_SEQ_NOT_FOUND')

        self.assertEqual(
            parse_master_dna_fastq.parse_dna_seq(seq_b_short, chr(50) * len(seq_b_short),
                chr(33), 7, 15),
            'END_SEQ_NOT_FOUND')

if __name__ == '__main__':
    unittest.main()


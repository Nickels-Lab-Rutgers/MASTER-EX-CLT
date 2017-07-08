#!/usr/bin/env python2.7
import env
import unittest
from src import rnautil

class TestSeqRmatch(unittest.TestCase):
    def test_seq_rmatch(self):
        self.assertTrue(rnautil.seq_rmatch('ACGT', 'ACGT'))
        self.assertTrue(rnautil.seq_rmatch('AACGT', 'ACGT'))
        self.assertTrue(rnautil.seq_rmatch('ACGT', 'AACGT'))
        self.assertTrue(rnautil.seq_rmatch('T', 'T'))
        self.assertTrue(rnautil.seq_rmatch('T', 'GT'))
        self.assertTrue(rnautil.seq_rmatch('GGGGGT', 'T'))

        self.assertFalse(rnautil.seq_rmatch('ACGA', 'ACGT'))
        self.assertFalse(rnautil.seq_rmatch('CGA', 'ACGT'))
        self.assertFalse(rnautil.seq_rmatch('TGT', 'ACGT'))
        self.assertFalse(rnautil.seq_rmatch('C', 'ACGT'))
        self.assertFalse(rnautil.seq_rmatch('CT', 'ACGT'))
        self.assertFalse(rnautil.seq_rmatch('CCGT', 'ACGT'))
        self.assertFalse(rnautil.seq_rmatch('AAGCGT', 'ACGT'))

    def test_check_match(self):
        self.assertTrue(rnautil.check_match('ACGT', 'ACGT'))
        self.assertTrue(rnautil.check_match('AACGT', 'ACGT'))
        self.assertTrue(rnautil.check_match('ACGT', 'AACGT'))
        self.assertTrue(rnautil.check_match('T', 'T'))
        self.assertTrue(rnautil.check_match('T', 'GT'))
        self.assertTrue(rnautil.check_match('GGGGGT', 'T'))

        self.assertFalse(rnautil.check_match('ACGA', 'ACGT'))
        self.assertFalse(rnautil.check_match('CGA', 'ACGT'))
        self.assertFalse(rnautil.check_match('TGT', 'ACGT'))
        self.assertFalse(rnautil.check_match('C', 'ACGT'))
        self.assertFalse(rnautil.check_match('CT', 'ACGT'))
        self.assertFalse(rnautil.check_match('CCGT', 'ACGT'))
        self.assertFalse(rnautil.check_match('AAGCGT', 'ACGT'))

if __name__ == '__main__':
    unittest.main()

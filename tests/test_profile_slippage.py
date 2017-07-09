#!/usr/bin/env python2.7
import env
import unittest
from src import profile_slippage as ps

class TestSlippageTypeCountTable(unittest.TestCase):
    def test_invalid_slippage_type(self):
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('ACGTACGG', 'H', 1)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AGGTTCGG', 'H', 2)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AGGTTCGG', 'H', 7)

        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('ACGTACGG', 'D', 6)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAGTACGG', 'D', 0)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('ACGTACGG', 'D', 7)

        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('ACGTACGG', 'I', 2)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAGTACGG', 'I', 0)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAGTACGG', 'I', 6)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAGTACGG', 'I', 7)

        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAGTACGG', 'k', 7)
        with self.assertRaises(ValueError):
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAGTACGG', '1', 7)

    def test_is_valid_slippage_type(self):
        # Invalid types
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('ACGTACGG', 'H', 1))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AGGTTCGG', 'H', 2))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AGGTTCGG', 'H', 7))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('ACGTACGG', 'D', 6))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGG', 'D', 0))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('ACGTACGG', 'D', 7))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('ACGTACGG', 'I', 2))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGG', 'I', 0))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGG', 'I', 6))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGC', 'I', 6))
        self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGG', 'I', 7))

        # Unimplemented types
        with self.assertRaises(ValueError):
            self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGG', 'o', 7))
        with self.assertRaises(ValueError):
            self.assertFalse(ps.SlippageTypeCountTable.is_valid_slippage_type('AAGTACGG', 'p', 7))

        # Valid types
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('CAATACGG', 'I', 0))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('CAAAACGG', 'I', 0))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('CAATACGG', 'I', 5))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('AAA', 'H', 0))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('CCC', 'H', 1))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('GGGGG', 'H', 3))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('CAATACGG', 'D', 0))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('CAATACGG', 'D', 2))
        self.assertTrue(ps.SlippageTypeCountTable.is_valid_slippage_type('AC', 'D', 0))

    def test_get_slippage_type_tup(self):
        self.assertEqual(ps.SlippageTypeCountTable('CAATACGG', 'I', 0).get_slippage_type_repr_tup(), 
                         ('CAATACGG', 'I', 0))
        self.assertEqual(ps.SlippageTypeCountTable('CAAAACGG', 'I', 0).get_slippage_type_repr_tup(), 
                         ('CAAAACGG', 'I', 0))
        self.assertEqual(ps.SlippageTypeCountTable('CAATACGG', 'I', 5).get_slippage_type_repr_tup(), 
                         ('CAATACGG', 'I', 5))
        self.assertEqual(ps.SlippageTypeCountTable('AAA', 'H', 0).get_slippage_type_repr_tup(), 
                         ('AAA', 'H', 0))
        self.assertEqual(ps.SlippageTypeCountTable('CCC', 'H', 1).get_slippage_type_repr_tup(), 
                         ('CCC', 'H', 1))
        self.assertEqual(ps.SlippageTypeCountTable('GGGGG', 'H', 3).get_slippage_type_repr_tup(), 
                         ('GGGGG', 'H', 3))
        self.assertEqual(ps.SlippageTypeCountTable('CAATACGG', 'D', 0).get_slippage_type_repr_tup(), 
                         ('CAATACGG', 'D', 0))
        self.assertEqual(ps.SlippageTypeCountTable('CAATACGG', 'D', 2).get_slippage_type_repr_tup(), 
                         ('CAATACGG', 'D', 2))
        self.assertEqual(ps.SlippageTypeCountTable('AC', 'D', 0).get_slippage_type_repr_tup(), 
                         ('AC', 'D', 0))

    
    def test_valid_slippage_type(self):
        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('CAATACGG', 'I', 0)
        except:
            self.fail("Valid slippage type failed to construct. ")

        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('CAAAACGG', 'I', 0)
        except:
            self.fail("Valid slippage type failed to construct. ")

        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('CAATACGG', 'I', 5)
        except:
            self.fail("Valid slippage type failed to construct. ")

        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AAA', 'H', 0)
        except:
            self.fail("Valid slippage type failed to construct. ")
        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('CCC', 'H', 1)
        except:
            self.fail("Valid slippage type failed to construct. ")
        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('GGGGG', 'H', 3)
        except:
            self.fail("Valid slippage type failed to construct. ")

        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('CAATACGG', 'D', 0)
        except:
            self.fail("Valid slippage type failed to construct. ")
        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('CAATACGG', 'D', 2)
        except:
            self.fail("Valid slippage type failed to construct. ")
        try:
            slp_type_cnt_tbl = ps.SlippageTypeCountTable('AC', 'D', 0)
        except:
            self.fail("Valid slippage type failed to construct. ")



    def test_is_slippage_seq(self):
        # Homopolymeric tranct
        h_cnt_tbl = ps.SlippageTypeCountTable('AACCTT', 'H', 0)
        self.assertTrue(h_cnt_tbl.is_slippage_seq('AAACCTT'))
        self.assertTrue(h_cnt_tbl.is_slippage_seq('AAAACCTT'))
        self.assertTrue(h_cnt_tbl.is_slippage_seq('AAAAACCTT'))
        self.assertTrue(h_cnt_tbl.is_slippage_seq('AAAAAACCTT'))
        self.assertTrue(h_cnt_tbl.is_slippage_seq('AAAAAAACCTT'))

        self.assertFalse(h_cnt_tbl.is_slippage_seq('ACCTT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('AACCTT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('AAACCTC'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('AAACCCT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('CAACCCT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('ATACCCT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('AAATACCCT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('CTT'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq('T'))
        self.assertFalse(h_cnt_tbl.is_slippage_seq(''))


        # Homopolymeric tranct internal
        i_cnt_tbl = ps.SlippageTypeCountTable('CAATACGG', 'I', 0)
        self.assertTrue(i_cnt_tbl.is_slippage_seq('CAAATACGG'))
        self.assertTrue(i_cnt_tbl.is_slippage_seq('CAAAATACGG'))
        self.assertTrue(i_cnt_tbl.is_slippage_seq('CAAAAATACGG'))

        self.assertFalse(i_cnt_tbl.is_slippage_seq('CAATACGG'))
        # mismatch at non-slippage region
        self.assertFalse(i_cnt_tbl.is_slippage_seq('CAAAATATGG'))
        self.assertFalse(i_cnt_tbl.is_slippage_seq('CAAATACGT'))
        # mismatch at slippage region
        self.assertFalse(i_cnt_tbl.is_slippage_seq('CTAATACGG'))
        # short
        self.assertFalse(i_cnt_tbl.is_slippage_seq('CAATACG'))
        self.assertFalse(i_cnt_tbl.is_slippage_seq('AATACG'))
        self.assertFalse(i_cnt_tbl.is_slippage_seq('ATACG'))
        self.assertFalse(i_cnt_tbl.is_slippage_seq('TACG'))
        self.assertFalse(i_cnt_tbl.is_slippage_seq(''))


        # Dinucleotide
        d_cnt_tbl = ps.SlippageTypeCountTable('ACGT', 'D', 0)
        self.assertTrue(d_cnt_tbl.is_slippage_seq('ACACGT'))
        self.assertTrue(d_cnt_tbl.is_slippage_seq('ACACACGT'))
        self.assertTrue(d_cnt_tbl.is_slippage_seq('ACACACACGT'))

        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACG'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('AC'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq(''))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACACGA'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('CACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('CACACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('CACACACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('CACACACACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACATACACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACAGATACGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACACCGT'))
        self.assertFalse(d_cnt_tbl.is_slippage_seq('ACACGCGT'))

    def test_h_slippage_add_rna_seq(self):
        # Homopolymeric tranct
        h_cnt_tbl = ps.SlippageTypeCountTable('AACCTT', 'H', 0)
        rna_seq1 = 'AAACCTT'
        h_cnt_tbl.add_rna_seq(rna_seq1, 2, 1)
        self.assertEqual(2, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(1, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        h_cnt_tbl.add_rna_seq(rna_seq1, 2, 1)
        self.assertEqual(4, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(2, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        h_cnt_tbl.add_rna_seq(rna_seq1, 2, 1)
        self.assertEqual(6, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(3, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])

        rna_seq2 = 'AAAACCTT'
        h_cnt_tbl.add_rna_seq(rna_seq2, 7, 4)
        self.assertEqual(7, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(4, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        h_cnt_tbl.add_rna_seq(rna_seq2, 7, 4)
        self.assertEqual(14, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(8, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        h_cnt_tbl.add_rna_seq(rna_seq2, 7, 4)
        self.assertEqual(21, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(12, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])

        rna_seq3 = 'AAAAACCTT'
        h_cnt_tbl.add_rna_seq(rna_seq3, 8, 5)
        self.assertEqual(8, h_cnt_tbl._slippage_cnt_dict[len(rna_seq3)]['raw_cnt'])
        self.assertEqual(5, h_cnt_tbl._slippage_cnt_dict[len(rna_seq3)]['tn_cnt'])
        h_cnt_tbl.add_rna_seq(rna_seq3, 8, 5)
        self.assertEqual(16, h_cnt_tbl._slippage_cnt_dict[len(rna_seq3)]['raw_cnt'])
        self.assertEqual(10, h_cnt_tbl._slippage_cnt_dict[len(rna_seq3)]['tn_cnt'])
        
        h_cnt_tbl.add_rna_seq('ACCTT', 20, 10)
        h_cnt_tbl.add_rna_seq('AACCTT', 20, 10)
        h_cnt_tbl.add_rna_seq('AAACCTC', 20, 10)
        h_cnt_tbl.add_rna_seq('AAACCCT', 20, 10)
        h_cnt_tbl.add_rna_seq('CAACCCT', 20, 10)
        h_cnt_tbl.add_rna_seq('ATACCCT', 20, 10)
        h_cnt_tbl.add_rna_seq('AAATACCT', 20, 10)
        h_cnt_tbl.add_rna_seq('AAATACCCT', 20, 10)
        h_cnt_tbl.add_rna_seq('CTT', 20, 10)
        h_cnt_tbl.add_rna_seq('T', 20, 10)
        h_cnt_tbl.add_rna_seq('', 20, 10)
        self.assertEqual(6, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(3, h_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        self.assertEqual(21, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(12, h_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        self.assertEqual(16, h_cnt_tbl._slippage_cnt_dict[len(rna_seq3)]['raw_cnt'])
        self.assertEqual(10, h_cnt_tbl._slippage_cnt_dict[len(rna_seq3)]['tn_cnt'])

        self.assertEqual(3, len(h_cnt_tbl._slippage_cnt_dict))

    def test_i_slippage_add_rna_seq(self):
        # Homopolymeric tranct internal
        i_cnt_tbl = ps.SlippageTypeCountTable('CAATACGG', 'I', 0)
        rna_seq1 = 'CAAATACGG'
        i_cnt_tbl.add_rna_seq(rna_seq1, 5, 3)
        self.assertEqual(5, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(3, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        i_cnt_tbl.add_rna_seq(rna_seq1, 5, 3)
        self.assertEqual(10, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(6, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        i_cnt_tbl.add_rna_seq(rna_seq1, 5, 3)
        self.assertEqual(15, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(9, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])

        rna_seq2 = 'CAAAATACGG'
        i_cnt_tbl.add_rna_seq(rna_seq2, 6, 2)
        self.assertEqual(6, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(2, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        i_cnt_tbl.add_rna_seq(rna_seq2, 6, 2)
        self.assertEqual(12, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(4, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        i_cnt_tbl.add_rna_seq(rna_seq2, 6, 2)
        self.assertEqual(18, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(6, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])

        # mismatch at non-slippage region
        i_cnt_tbl.add_rna_seq('CAAAATATGG', 20, 10)
        i_cnt_tbl.add_rna_seq('CAAATACGT', 20, 10)
        # mismatch at slippage region
        i_cnt_tbl.add_rna_seq('CTAATACGG', 20, 10)
        # short
        i_cnt_tbl.add_rna_seq('CAATACG', 20, 10)
        i_cnt_tbl.add_rna_seq('AATACG', 20, 10)
        i_cnt_tbl.add_rna_seq('ATACG', 20, 10)
        i_cnt_tbl.add_rna_seq('TACG', 20, 10)
        i_cnt_tbl.add_rna_seq('', 20, 10)

        self.assertEqual(15, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(9, i_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        self.assertEqual(18, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(6, i_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])

        self.assertEqual(2, len(i_cnt_tbl._slippage_cnt_dict))

    def test_d_slippage_add_rna_seq(self):
        # Dinucleotide
        d_cnt_tbl = ps.SlippageTypeCountTable('ACGT', 'D', 0)
        rna_seq1 = 'ACACGT'
        d_cnt_tbl.add_rna_seq(rna_seq1, 3, 1)
        self.assertEqual(3, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(1, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        d_cnt_tbl.add_rna_seq(rna_seq1, 3, 1)
        self.assertEqual(6, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(2, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        d_cnt_tbl.add_rna_seq(rna_seq1, 3, 1)
        self.assertEqual(9, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(3, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])

        rna_seq2 = 'ACACACGT'
        d_cnt_tbl.add_rna_seq(rna_seq2, 5, 2)
        d_cnt_tbl.add_rna_seq(rna_seq2, 5, 2)
        self.assertEqual(10, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(4, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        d_cnt_tbl.add_rna_seq(rna_seq2, 3, 1)
        self.assertEqual(13, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(5, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])
        d_cnt_tbl.add_rna_seq(rna_seq2, 7, 5)
        self.assertEqual(20, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(10, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])

        d_cnt_tbl.add_rna_seq('ACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('ACG', 20, 10)
        d_cnt_tbl.add_rna_seq('AC', 20, 10)
        d_cnt_tbl.add_rna_seq('', 20, 10)
        d_cnt_tbl.add_rna_seq('ACACGA', 20, 10)
        d_cnt_tbl.add_rna_seq('CACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('CACACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('CACACACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('CACACACACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('ACATACACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('ACAGATACGT', 20, 10)
        d_cnt_tbl.add_rna_seq('ACACCGT', 20, 10)
        d_cnt_tbl.add_rna_seq('ACACGCGT', 20, 10)

        self.assertEqual(9, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['raw_cnt'])
        self.assertEqual(3, d_cnt_tbl._slippage_cnt_dict[len(rna_seq1)]['tn_cnt'])
        self.assertEqual(20, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['raw_cnt'])
        self.assertEqual(10, d_cnt_tbl._slippage_cnt_dict[len(rna_seq2)]['tn_cnt'])

        self.assertEqual(2, len(d_cnt_tbl._slippage_cnt_dict))


class TestDNASlippageTable(unittest.TestCase):
    def test_dna_sequence_check(self):
        with self.assertRaises(ValueError):
            ps.DNASlippageTable('')
        with self.assertRaises(ValueError):
            ps.DNASlippageTable('A')
        try:
            ps.DNASlippageTable('AA')
        except:
            self.fail('dna_prmt_seq with two nt should be ok')

    def test_get_slippage_type_tuple(self):
        def test_slippage_type_tuple_equal(x, y):
            set_x = set(x)
            set_y = set(y)
            self.assertEqual(len(x), len(set_x))
            self.assertEqual(len(y), len(set_y))
            self.assertEqual(set_x, set_y)

        test_slippage_type_tuple_equal(ps.DNASlippageTable.get_slippage_type_tuple('AA'),
                                       (('H', 0), ))
        test_slippage_type_tuple_equal(ps.DNASlippageTable.get_slippage_type_tuple('AC'),
                                       (('D', 0), ))
        test_slippage_type_tuple_equal(ps.DNASlippageTable.get_slippage_type_tuple('ACC'),
                                       (('D', 0), ('I', 0), ('H', 1)))
        test_slippage_type_tuple_equal(ps.DNASlippageTable.get_slippage_type_tuple('ACCC'),
                                       (('D', 0), ('I', 0), ('H', 1), ('H', 2)))
        test_slippage_type_tuple_equal(ps.DNASlippageTable.get_slippage_type_tuple('ACCCG'),
                                       (('D', 0), ('I', 0), ('H', 1), ('H', 2), ('D', 3)))
        test_slippage_type_tuple_equal(ps.DNASlippageTable.get_slippage_type_tuple('ACCGCC'),
                                       (('D', 0), ('D', 2), ('D', 3), ('I', 0), ('I', 3), ('H', 1), ('H', 4)))

    def test_slp_type_cnt_tbl_tup(self):
        def assert_slp_type_cnt_tbl_tup_equal(dna_prmt_seq, slp_type_tup):
            dna_slp_tbl = ps.DNASlippageTable(dna_prmt_seq)
            dna_slp_tbl_slp_typ_repr_tup = tuple(map(lambda st: st.get_slippage_type_repr_tup(), 
                                                     dna_slp_tbl._slp_type_cnt_tbl_dict.values()))
            slp_type_cnt_tbl_tup_repr_lst = map(lambda t: (dna_prmt_seq,) + t, slp_type_tup)
            self.assertEqual(len(dna_slp_tbl_slp_typ_repr_tup), len(set(dna_slp_tbl_slp_typ_repr_tup)))
            self.assertEqual(set(dna_slp_tbl_slp_typ_repr_tup), set(slp_type_cnt_tbl_tup_repr_lst))

        assert_slp_type_cnt_tbl_tup_equal('AA', (('H', 0), ))
        assert_slp_type_cnt_tbl_tup_equal('AC', 
                                          (('D', 0), ))
        assert_slp_type_cnt_tbl_tup_equal('ACC', 
                                          (('D', 0), ('I', 0), ('H', 1)))
        assert_slp_type_cnt_tbl_tup_equal('ACCC', 
                                          (('D', 0), ('I', 0), ('H', 1), ('H', 2)))
        assert_slp_type_cnt_tbl_tup_equal('ACCCG', 
                                          (('D', 0), ('I', 0), ('H', 1), ('H', 2), ('D', 3)))
        assert_slp_type_cnt_tbl_tup_equal('ACCGCC', 
                                          (('D', 0), ('D', 2), ('D', 3), ('I', 0), ('I', 3), ('H', 1), ('H', 4)))

    def test_add_rna_seq(self):
        dna_slp_tbl = ps.DNASlippageTable('ACCGCC')
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 0)
        with self.assertRaises(ValueError):
            dna_slp_tbl.add_rna_seq('', 1, 1)
        with self.assertRaises(ValueError):
            dna_slp_tbl.add_rna_seq('AACCGCC', 1, 1)
        with self.assertRaises(ValueError):
            dna_slp_tbl.add_rna_seq('AAACCGCC', 1, 1)

        rna_seq1 = 'ACC'
        dna_slp_tbl.add_rna_seq(rna_seq1, 10, 3)
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 1)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 3)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 10)
        for slp_type_cnt_tbl in dna_slp_tbl._slp_type_cnt_tbl_dict.values():
            self.assertEqual(len(slp_type_cnt_tbl._slippage_cnt_dict), 0)

        dna_slp_tbl.add_rna_seq(rna_seq1, 9, 2)
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 1)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 5)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 19)
        for slp_type_cnt_tbl in dna_slp_tbl._slp_type_cnt_tbl_dict.values():
            self.assertEqual(len(slp_type_cnt_tbl._slippage_cnt_dict), 0)

        rna_seq2 = 'GCC'
        assert len(rna_seq1) == len(rna_seq2)
        dna_slp_tbl.add_rna_seq(rna_seq2, 10, 3)
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 1)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 5)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 19)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq2)]['matched_tn_cnt'], 3)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq2)]['matched_raw_cnt'], 10)

        dna_slp_tbl.add_rna_seq(rna_seq2, 10, 3)
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 1)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 5)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 19)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq2)]['matched_tn_cnt'], 6)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq2)]['matched_raw_cnt'], 20)

        for slp_type_cnt_tbl in dna_slp_tbl._slp_type_cnt_tbl_dict.values():
            self.assertEqual(len(slp_type_cnt_tbl._slippage_cnt_dict), 0)

        rna_seq3 = 'CCC'
        assert len(rna_seq3) == len(rna_seq1)
        dna_slp_tbl.add_rna_seq(rna_seq3, 3, 2)
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 1)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 7)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 22)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['matched_tn_cnt'], 6)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['matched_raw_cnt'], 20)

        for slp_type_tup in dna_slp_tbl._slp_type_cnt_tbl_dict:
            if slp_type_tup != ('H', 4):
                self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict), 0)
            self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[('H', 4)]._slippage_cnt_dict), 1)
            self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[('H', 4)]._slippage_cnt_dict[len(rna_seq3)]['tn_cnt'], 2)
            self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[('H', 4)]._slippage_cnt_dict[len(rna_seq3)]['raw_cnt'], 3)

        dna_slp_tbl.add_rna_seq(rna_seq3, 3, 2)
        self.assertEqual(len(dna_slp_tbl._rna_cnt_dict), 1)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 9)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 25)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['matched_tn_cnt'], 6)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['matched_raw_cnt'], 20)

        for slp_type_tup in dna_slp_tbl._slp_type_cnt_tbl_dict:
            if slp_type_tup != ('H', 4):
                self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict), 0)
            self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[('H', 4)]._slippage_cnt_dict), 1)
            self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[('H', 4)]._slippage_cnt_dict[len(rna_seq3)]['tn_cnt'], 4)
            self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[('H', 4)]._slippage_cnt_dict[len(rna_seq3)]['raw_cnt'], 6)

        rna_seq4 = 'GCCC'
        dna_slp_tbl.add_rna_seq(rna_seq4, 5, 3)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_tn_cnt'], 9)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['mismatched_raw_cnt'], 25)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['matched_tn_cnt'], 6)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq1)]['matched_raw_cnt'], 20)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq4)]['mismatched_tn_cnt'], 3)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq4)]['mismatched_raw_cnt'], 5)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq4)]['matched_tn_cnt'], 0)
        self.assertEqual(dna_slp_tbl._rna_cnt_dict[len(rna_seq4)]['matched_raw_cnt'], 0)
        
        for slp_type_tup in dna_slp_tbl._slp_type_cnt_tbl_dict:
            if slp_type_tup == ('H', 4):
                self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict), 1)
                self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict[len(rna_seq3)]['tn_cnt'], 4)
                self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict[len(rna_seq3)]['raw_cnt'], 6)
            elif slp_type_tup == ('I', 3):
                self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict), 1)
                self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict[len(rna_seq4)]['tn_cnt'], 3)
                self.assertEqual(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict[len(rna_seq4)]['raw_cnt'], 5)
            else:
                self.assertEqual(len(dna_slp_tbl._slp_type_cnt_tbl_dict[slp_type_tup]._slippage_cnt_dict), 0)
        
if __name__ == '__main__':
    unittest.main()
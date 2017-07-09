#!/usr/bin/env python2.7
import argparse
import rnautil
import collections

class DNASlippageTable(object):
    """docstring for DNASlippageTable"""

    def __init__(self, dna_prmt_seq):
        super(DNASlippageTable, self).__init__()
        if len(dna_prmt_seq) < 2:
            msg = 'DNA promoter sequence should have at least 2 nucleotides: ' + dna_prmt_seq
            raise ValueError(msg)
        
        self._dna_prmt_seq = dna_prmt_seq
        dna_ps_slp_type_tup_tup = self.get_slippage_type_tuple(dna_prmt_seq)
        self._slp_type_cnt_tbl_dict = dict(zip(dna_ps_slp_type_tup_tup,
                                               tuple(SlippageTypeCountTable(dna_prmt_seq, 
                                                                            slp_type_tup[0], 
                                                                            slp_type_tup[1]) 
                                                     for slp_type_tup in dna_ps_slp_type_tup_tup)))

        # {length : {'matched_tn_cnt' : 0, ...}}
        self._rna_cnt_dict = collections.defaultdict(lambda : {'matched_tn_cnt' : 0,
                                                               'matched_raw_cnt' : 0,
                                                               'mismatched_tn_cnt' : 0,
                                                               'mismatched_raw_cnt' : 0})

    @staticmethod
    def get_slippage_type_tuple(dna_prmt_seq):
        slippage_type_list = []
        # Implemented: H, D, and I
        # I match starts at the cap. Z in ZAA is cap.
        for i in xrange(len(dna_prmt_seq) - 1):
            if dna_prmt_seq[i] == dna_prmt_seq[i + 1]:
                slippage_type_list.append(('H', i))
                if i - 1 >= 0 and dna_prmt_seq[i - 1] != dna_prmt_seq[i]:
                    slippage_type_list.append(('I', i - 1))
            else:
                slippage_type_list.append(('D', i))
        return tuple(slippage_type_list)

    def add_rna_seq(self, rna_prmt_seq, raw_cnt, tn_cnt):
        if len(rna_prmt_seq) > len(self._dna_prmt_seq):
            raise ValueError('rna_prmt_seq {} is \
longer than dna_prmt_seq {}.'.format(rna_prmt_seq, self._dna_prmt_seq))
        
        if len(rna_prmt_seq) == 0:
            raise ValueError('rna_prmt_seq is empty. DNA: {}. \
RNA: {}. Raw count: {}. TN count: {}'.format(self._dna_prmt_seq, rna_prmt_seq, raw_cnt, tn_cnt))
        
        if rnautil.seq_rmatch(self._dna_prmt_seq, rna_prmt_seq):
            self._rna_cnt_dict[len(rna_prmt_seq)]['matched_tn_cnt'] += tn_cnt
            self._rna_cnt_dict[len(rna_prmt_seq)]['matched_raw_cnt'] += raw_cnt
        else:
            self._rna_cnt_dict[len(rna_prmt_seq)]['mismatched_tn_cnt'] += tn_cnt
            self._rna_cnt_dict[len(rna_prmt_seq)]['mismatched_raw_cnt'] += raw_cnt

        for slp_cnt_tbl in self._slp_type_cnt_tbl_dict.values():
            slp_cnt_tbl.add_rna_seq(rna_prmt_seq, raw_cnt, tn_cnt)


# Count table for each slippage type. 
# This class only records the count of slippage reads. 
# General counts are stored in @class DNASlippageTable. 
# 
# Implemented patterns:
# 1. H: homopolymeric tract, NAAN. At least two consecutive same nt started at start ind. 
#       Slippage reads have extra homopolymeric sequence upstream of start ind. A...AAN
# 2. I: homopolymeric tract internal, ZAAN. Z means different from A. 
#       At least two consecutive same nt started at start ind + 1, 
#       and a different nt at ind. Slippage reads have extra homopolymeric sequence upstream of
#       ind + 1 that is capped by Z. 
# 3. D: dinucleotide, NACN. Two different nts start at ind. Slippage reads have extra dinucleotide
#       repeats upstream of ind. 
class SlippageTypeCountTable(object):
    """docstring for SlippageTypeCountTable"""

    def __init__(self, dna_prmt_seq, slippage_type, match_start_ind):
        super(SlippageTypeCountTable, self).__init__()
        if not self.is_valid_slippage_type(dna_prmt_seq, slippage_type, match_start_ind):
            msg = "Invalid slippage type: " + str((dna_prmt_seq, slippage_type, match_start_ind))
            raise ValueError(msg)
        
        self._dna_prmt_seq = dna_prmt_seq
        self._slippage_type = slippage_type
        self._match_start_ind = match_start_ind
        self._matched_rna_seq = dna_prmt_seq[match_start_ind:]
        # {} is 'constructor', so it is ok to use this lambda as factory function.
        # {length : {'tn_cnt' : 0, 'raw_cnt' : 0}, ...}
        self._slippage_cnt_dict = collections.defaultdict(lambda : {'tn_cnt' : 0, 'raw_cnt' : 0})

    def get_slippage_type_repr_tup(self):
        return (self._dna_prmt_seq, self._slippage_type, self._match_start_ind)

    def add_rna_seq(self, rna_prmt_seq, raw_cnt, tn_cnt):
        if self.is_slippage_seq(rna_prmt_seq):
            self._slippage_cnt_dict[len(rna_prmt_seq)]['tn_cnt'] += tn_cnt
            self._slippage_cnt_dict[len(rna_prmt_seq)]['raw_cnt'] += raw_cnt
    
    @staticmethod
    def is_valid_slippage_type(dna_prmt_seq, slippage_type, match_start_ind):
        if slippage_type == 'H':
            if (match_start_ind + 1 < len(dna_prmt_seq)
                and dna_prmt_seq[match_start_ind] == dna_prmt_seq[match_start_ind + 1]):
                return True
        elif slippage_type == 'I':
            if (match_start_ind + 2 < len(dna_prmt_seq)
                and dna_prmt_seq[match_start_ind] != dna_prmt_seq[match_start_ind + 1]
                and dna_prmt_seq[match_start_ind + 1] == dna_prmt_seq[match_start_ind + 2]):
                return True
        elif slippage_type == 'D':
            if (match_start_ind + 1 < len(dna_prmt_seq)
                and dna_prmt_seq[match_start_ind] != dna_prmt_seq[match_start_ind + 1]):
                return True
        else:
            raise ValueError('Unimplemented slippage type: ' + slippage_type)
        return False

    def is_slippage_seq(self, rna_prmt_seq):
        # Sequence shorter than matched seq is not slippage
        if len(rna_prmt_seq) <= len(self._matched_rna_seq):
            return False

        # Example:
        # Homopolymeric tract. Matching starts at 3. 
        #                  DNA: CGTAAAA
        #          match start:    ^
        #                  RNA:  AAAAAA
        # rna_non_slippage_seq:    ****
        #     rna_slippage_seq:  **
        #
        # Homopolymeric tract internal. Matching starts at 2.
        #                  DNA: CGTAAAA
        #          match start:   ^
        #                  RNA: TAAAAAA
        # rna_non_slippage_seq:    ****
        #     rna_slippage_seq: ***
        #
        # Dinucleotide. Matching starts at 2. 
        #                  DNA: CGTAAAA
        #          match start:   ^
        #                  RNA: TATAAAA
        # rna_non_slippage_seq:   *****
        #     rna_slippage_seq: **
        if self._slippage_type == 'H':
            rna_non_slippage_seq_len = len(self._matched_rna_seq)
            rna_non_slippage_seq = rna_prmt_seq[-rna_non_slippage_seq_len:]
            rna_slippage_seq = rna_prmt_seq[:-rna_non_slippage_seq_len]

            if rna_non_slippage_seq != self._matched_rna_seq:
                return False
            
            # homopolymeric tract nucleotide
            ht_nt = self._dna_prmt_seq[self._match_start_ind]
            for rna_nt in rna_slippage_seq:
                if rna_nt != ht_nt:
                    return False
        elif self._slippage_type == 'I':
            rna_non_slippage_seq_len = len(self._matched_rna_seq) - 1
            rna_non_slippage_seq = rna_prmt_seq[-rna_non_slippage_seq_len:]
            rna_slippage_seq = rna_prmt_seq[:-rna_non_slippage_seq_len]

            if rna_non_slippage_seq != self._matched_rna_seq[1:]:
                return False

            # homopolymeric tract nucleotide
            ht_nt = self._dna_prmt_seq[self._match_start_ind + 1]
            cap_nt = self._dna_prmt_seq[self._match_start_ind]
            if rna_slippage_seq[0] != cap_nt:
                return False

            for rna_nt in rna_slippage_seq[1:]:
                if rna_nt != ht_nt:
                    return False
        elif self._slippage_type == 'D':
            rna_non_slippage_seq_len = len(self._matched_rna_seq)
            rna_non_slippage_seq = rna_prmt_seq[-rna_non_slippage_seq_len:]
            rna_slippage_seq = rna_prmt_seq[:-rna_non_slippage_seq_len]

            if rna_non_slippage_seq != self._matched_rna_seq:
                return False

            if len(rna_slippage_seq) % 2 != 0:
                return False

            dinuc_seq = self._dna_prmt_seq[self._match_start_ind:self._match_start_ind + 2]
            for rna_dinuc_ind in xrange(0, len(rna_slippage_seq), 2):
                if rna_slippage_seq[rna_dinuc_ind:rna_dinuc_ind + 2] != dinuc_seq:
                    return False
        else:
            msg = 'Violation of invariant: unimplemented slippage type {}'.format(self._slippage_type)
            raise ValueError(msg)

        return True


# RNA length ranges fraom 1 to (dna_len + 5). Inclusive.
class LibrarySlippageTable(object):
    """docstring for LibrarySlippageTable"""
    _MAX_NUM_BP_UPSTREM_RND_REGION = 5

    def __init__(self, rna_parsed_tbl_iterator, beg_seq):
        super(LibrarySlippageTable, self).__init__()

        if len(beg_seq) < self._MAX_NUM_BP_UPSTREM_RND_REGION:
            msg = 'Beginning sequence should have at \
least {} nucleotides: {}'.format(self._MAX_NUM_BP_UPSTREM_RND_REGION, beg_seq)
            raise ValueError(msg)

        self._beg_seq = beg_seq
        self.slippage_cnt_dict = self.get_slippage_cnt_dict(rna_parsed_tbl_iterator)

    def get_slippage_cnt_dict(self, rna_parsed_tbl_iterator):
        # {dna : slippage_counter}

        slpg_cnt_dict = {}
        #for rna_parsed_record in rna_parsed_tbl_iterator:
        #    complete_dna_prmt_seq = rna_parsed_record.dna_prmt_seq
        pass



    def get_ptn_slippage_tbl(self, ptn_list):
        pass


class PtnSlippageTbl(object):
    """docstring for PtnSlippageTbl"""
    def __init__(self):
        super(PtnSlippageTbl, self).__init__()

    def write_slippage_tbl(self, output_fn):
        pass
        

def nt_sequence(nt_seq):
    if rnautil.is_valid_seq(nt_seq):
        return nt_seq
    else:
        msg = "Invalid nucleotide sequence: " + nt_seq
        raise argparse.ArgumentTypeError(msg)

def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('beg_seq', metavar = '<Beginning fixed sequence>', 
                            type = nt_sequence)

    arg_parser.add_argument('ptn_fn', metavar = '<slippage pattern table>', 
                            type = argparse.FileType('r'))

    arg_parser.add_argument('rna_parsed_tbl_iterator', 
                            metavar = '<rna parsed file>', 
                            type = rnautil.iterate_rna_parsed_record,
                            help = 'All DNA/RNA promoter region sequence count file. '
                                   'Tab delimited. '
                                   'Columns from left to right are RNA, DNA, XX, raw count, tn count, and XX. '
                                   'XX denotes deprecated fields. ')

    arg_parser.add_argument('ofile', metavar = '<output file>', type = argparse.FileType('w'))

    arg_parser.parse_args()

if __name__ == '__main__':
    main()
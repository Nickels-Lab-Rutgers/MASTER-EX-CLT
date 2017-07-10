#!/usr/bin/env python2.7
import argparse
import rnautil
import collections

# Fields started with '_' are read only. 

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

        self._beg_seq = beg_seq[-self._MAX_NUM_BP_UPSTREM_RND_REGION:]
        self._lib_slpg_cnt_dict = self.get_lib_slpg_cnt_dict(rna_parsed_tbl_iterator)
        # DNA promoter region sequences are required to have the same length. 
        self._cmpl_dna_prmt_seq_len = len(self._lib_slpg_cnt_dict.keys()[0])
        
    def get_lib_slpg_cnt_dict(self, rna_parsed_tbl_iterator):
        # {dna : dna_slippage_count_table}
        lib_slpg_cnt_dict = {}
        for rna_parsed_record in rna_parsed_tbl_iterator:
            # Only process rna_prmt_seq with less than _MAX_NUM_BP_UPSTREM_RND_REGION upstream
            # of dna_prmt_seq
            if (len(rna_parsed_record.rna_prmt_seq) 
                > self._MAX_NUM_BP_UPSTREM_RND_REGION + len(rna_parsed_record.dna_prmt_seq)):
                continue

            cmpl_dna_prmt_seq = self._beg_seq + rna_parsed_record.dna_prmt_seq 
            if cmpl_dna_prmt_seq not in lib_slpg_cnt_dict:
                lib_slpg_cnt_dict[cmpl_dna_prmt_seq] = DNASlippageTable(cmpl_dna_prmt_seq)

            lib_slpg_cnt_dict[cmpl_dna_prmt_seq].add_rna_seq(rna_parsed_record.rna_prmt_seq, 
                                                             rna_parsed_record.raw_cnt, 
                                                             rna_parsed_record.tn_cnt)
        return lib_slpg_cnt_dict

    def get_ptn_slippage_tbl(self, ptn_tup_tup):
        # all pattern slippage table
        a_ptn_slippage_tbl = PtnSlippageTbl()

        for ptn_tup in ptn_tup_tup:
            # single pattern cnt dict
            s_ptn_cnt_dict = dict(zip(xrange(1, self._cmpl_dna_prmt_seq_len + 1), 
                                      ({'matched_tn_cnt' : 0, 'matched_raw_cnt' : 0,
                                        'mismatched_tn_cnt' : 0, 'mismatched_raw_cnt' : 0,
                                        'slippage_tn_cnt' : 0, 'slippage_raw_cnt' : 0}
                                        for i in xrange(1, self._cmpl_dna_prmt_seq_len + 1))))
            num_matched_dna_prmt_seq = 0

            for cmpl_dna_prmt_seq in self._lib_slpg_cnt_dict:
                dna_prmt_seq = cmpl_dna_prmt_seq[self._MAX_NUM_BP_UPSTREM_RND_REGION:]
                if self.ptn_dna_match(ptn_tup, dna_prmt_seq):
                    num_matched_dna_prmt_seq += 1
                    dna_slippage_cnt_tbl = self._lib_slpg_cnt_dict[cmpl_dna_prmt_seq]
                    for length in dna_slippage_cnt_tbl._rna_cnt_dict:
                        s_ptn_cnt_dict[length]['matched_tn_cnt'] \
                            += dna_slippage_cnt_tbl._rna_cnt_dict[length]['matched_tn_cnt']
                        s_ptn_cnt_dict[length]['matched_raw_cnt'] \
                            += dna_slippage_cnt_tbl._rna_cnt_dict[length]['matched_raw_cnt']
                        s_ptn_cnt_dict[length]['mismatched_tn_cnt'] \
                            += dna_slippage_cnt_tbl._rna_cnt_dict[length]['mismatched_tn_cnt']
                        s_ptn_cnt_dict[length]['mismatched_raw_cnt'] \
                            += dna_slippage_cnt_tbl._rna_cnt_dict[length]['mismatched_raw_cnt']

                    if (ptn_tup.slippage_type, ptn_tup.match_start_ind + self._MAX_NUM_BP_UPSTREM_RND_REGION) \
                        not in dna_slippage_cnt_tbl._slp_type_cnt_tbl_dict:
                        raise ValueError('Slippage pattern not in dna_cnt_tbl: {}'.format(ptn_tup))

                    dna_slpg_type_cnt_tbl = dna_slippage_cnt_tbl._slp_type_cnt_tbl_dict[(ptn_tup.slippage_type, 
                                                    ptn_tup.match_start_ind + self._MAX_NUM_BP_UPSTREM_RND_REGION)]
                    for length in dna_slpg_type_cnt_tbl._slippage_cnt_dict:
                        s_ptn_cnt_dict[length]['slippage_tn_cnt'] \
                            += dna_slpg_type_cnt_tbl._slippage_cnt_dict[length]['tn_cnt']
                        s_ptn_cnt_dict[length]['slippage_raw_cnt'] \
                            += dna_slpg_type_cnt_tbl._slippage_cnt_dict[length]['raw_cnt']

            a_ptn_slippage_tbl.insert_s_ptn_slp_cnt_dict(ptn_tup, s_ptn_cnt_dict, num_matched_dna_prmt_seq)

        return a_ptn_slippage_tbl

    # Assume that slp_ptn_tup is valid.
    # pattern only has ACGTNZ
    # valid length
    @staticmethod
    def ptn_dna_match(slp_ptn_tup, dna_prmt_seq):
        ptn_seq = slp_ptn_tup.ptn_seq
        slp_type = slp_ptn_tup.slippage_type
        match_start_ind = slp_ptn_tup.match_start_ind
        
        if len(ptn_seq) != len(dna_prmt_seq):
            msg = "Length of pattern '{}' is different from length of DNA '{}'".format(ptn_seq, dna_prmt_seq)
            raise ValueError(msg)

        if slp_type == 'H':
            for i in xrange(len(ptn_seq)):
                if ptn_seq[i] != 'N':
                    if ptn_seq[i] == 'Z':
                        if dna_prmt_seq[i] == ptn_seq[match_start_ind]:
                            return False
                    elif ptn_seq[i] in 'ACGT':
                        if dna_prmt_seq[i] != ptn_seq[i]:
                            return False
                    else:
                        raise ValueError('Unimplemented pattern character: {}'.format(ptn_seq[i]))
        elif slp_type == 'I':
            for i in xrange(len(ptn_seq)):
                if ptn_seq[i] != 'N':
                    if ptn_seq[i] == 'Z':
                        if dna_prmt_seq[i] == ptn_seq[match_start_ind + 1]:
                            return False
                    elif ptn_seq[i] in 'ACGT':
                        if dna_prmt_seq[i] != ptn_seq[i]:
                            return False
                    else:
                        raise ValueError('Unimplemented pattern character: {}'.format(ptn_seq[i]))
        elif slp_type == 'D':
            for i in xrange(len(ptn_seq)):
                if ptn_seq[i] != 'N':
                    if ptn_seq[i] == 'Z':
                        raise ValueError('Invalid pattern: {}'.format(slp_ptn_tup))
                    elif ptn_seq[i] in 'ACGT':
                        if dna_prmt_seq[i] != ptn_seq[i]:
                            return False
                    else:
                        raise ValueError('Unimplemented pattern character: {}'.format(ptn_seq[i]))
        else:
            raise ValueError('Unimplemented slippage type: {}. {}'.format(slp_type, slp_ptn_tup))

        return True


# This class only stores the data and format the output.
# It also includes helper function to validate pattern. 
# Pattern counting is done by LibrarySlippageTable.
class PtnSlippageTbl(object):
    """docstring for PtnSlippageTbl"""
    def __init__(self):
        super(PtnSlippageTbl, self).__init__()
        # {ptn_tup : {length : {'tnXXcnt' : 0, 'rawXXcnt' : 0, ...}, ...}, ...}
        self._ptn_slp_cnt_dict = {}
        self._ptn_slp_cnt_dict_ordered_key_list = []
        self._num_dna_prmt_list = []

    def insert_s_ptn_slp_cnt_dict(self, ptn_tup, s_ptn_slp_cnt_dict, num_dna_prmt):
        if ptn_tup in self._ptn_slp_cnt_dict:
            raise ValueError('{} is duplicated.'.format(ptn_tup))
        
        self._ptn_slp_cnt_dict[ptn_tup] = s_ptn_slp_cnt_dict
        self._ptn_slp_cnt_dict_ordered_key_list.append(ptn_tup)
        self._num_dna_prmt_list.append(num_dna_prmt)
    
    # This requires that the length range of each slippage pattern to be the same.
    def fmt_slippage_tbl(self):
        slp_tbl_str = self.get_table_header()
        for ind, ptn_tup in enumerate(self._ptn_slp_cnt_dict_ordered_key_list):
            slp_tbl_str += self.fmt_ptn_tup(ptn_tup) + '\t'
            slp_tbl_str += str(self._num_dna_prmt_list[ind]) + '\t'
            s_ptn_slp_cnt_dict = self._ptn_slp_cnt_dict[ptn_tup]
            slp_tbl_str += self.fmt_s_ptn_slp_cnt_dict(s_ptn_slp_cnt_dict)
        return slp_tbl_str

    def has_same_len_range(self):
        a_ptn_len_list = map(lambda d: sorted(d.keys()), self._ptn_slp_cnt_dict.values())
        if len(a_ptn_len_list) == 0:
            return True

        for s_ptn_len_list in a_ptn_len_list[1:]:
            if s_ptn_len_list != a_ptn_len_list[0]:
                return False
        return True

    def get_table_header(self):
        if not self.has_same_len_range():
            raise ValueError('Pattern slippage table has different lengths.')

        if len(self._ptn_slp_cnt_dict) == 0:
            return ''

        len_list = sorted(self._ptn_slp_cnt_dict.values()[0].keys(), reverse = True)
        header_str = 'Pattern\tnumDNAPrmt\t'
        header_str += '\t'.join(map(lambda l: 'tnSlp{}bp'.format(l), len_list)) + '\t'
        header_str += '\t'.join(map(lambda l: 'rawSlp{}bp'.format(l), len_list)) + '\t'
        header_str += '\t'.join(map(lambda l: 'tnM{}bp'.format(l), len_list)) + '\t'
        header_str += '\t'.join(map(lambda l: 'rawM{}bp'.format(l), len_list)) + '\t'
        header_str += '\t'.join(map(lambda l: 'tnA{}bp'.format(l), len_list)) + '\t'
        header_str += '\t'.join(map(lambda l: 'rawA{}bp'.format(l), len_list)) + '\n'
        return header_str

    # From left to right
    # length from long to short
    # TNSLP, SLP, TNM, M, TNA, A
    @staticmethod
    def fmt_s_ptn_slp_cnt_dict(s_ptn_slp_cnt_dict):
        s_ptn_slp_cnt_tup = tuple(s_ptn_slp_cnt_dict[length] 
                                  for length in sorted(s_ptn_slp_cnt_dict.keys(), reverse = True))
        tn_slp_list = map(lambda cnt_dict: cnt_dict['slippage_tn_cnt'], 
                          s_ptn_slp_cnt_tup)
        raw_slp_list = map(lambda cnt_dict: cnt_dict['slippage_raw_cnt'], 
                           s_ptn_slp_cnt_tup)
        tn_m_list = map(lambda cnt_dict: cnt_dict['matched_tn_cnt'], 
                        s_ptn_slp_cnt_tup)
        raw_m_list = map(lambda cnt_dict: cnt_dict['matched_raw_cnt'], 
                         s_ptn_slp_cnt_tup)
        tn_a_list = map(lambda cnt_dict: cnt_dict['mismatched_tn_cnt'] + cnt_dict['matched_tn_cnt'], 
                        s_ptn_slp_cnt_tup)
        raw_a_list = map(lambda cnt_dict: cnt_dict['mismatched_raw_cnt'] + cnt_dict['matched_raw_cnt'], 
                         s_ptn_slp_cnt_tup)

        s_ptn_slp_cnt_str = ''
        for l in (tn_slp_list, raw_slp_list, tn_m_list, raw_m_list, tn_a_list, raw_a_list):
            s_ptn_slp_cnt_str += '\t'.join(map(str, l)) + '\t'

        s_ptn_slp_cnt_str += '\n'
        return s_ptn_slp_cnt_str

    
    @staticmethod
    def fmt_ptn_tup(slp_ptn_tup):
        return slp_ptn_tup.ptn_seq + '-' + slp_ptn_tup.slippage_type + str(slp_ptn_tup.match_start_ind)

    # slp_ptn is a tuple (ptn, slp_type, match_start_ind) that indicates the 
    # pattern of slippage. Examples,
    # ('NNAAANN', 'H', 2)
    # ('NZAAANN', 'H', 2)
    # ('NZAAANN', 'I', 1)
    # ('NTAAANN', 'D', 1)
    # N: [ACGT]
    # Z: nt different from the homopolymeric tract
    @staticmethod
    def is_valid_slp_ptn_tup(slp_ptn_tup):
        ptn_seq = slp_ptn_tup.ptn_seq
        slp_type = slp_ptn_tup.slippage_type
        match_start_ind = slp_ptn_tup.match_start_ind

        for ptn_nt in ptn_seq:
            if ptn_nt not in 'ACGTNZ':
                raise ValueError('Unimplemented pattern character: {}'.format(ptn_nt))
        
        if slp_type == 'H':
            # Implicit that ptn_seq[match_start_ind + 1] not in 'NZ'
            if (match_start_ind + 1 < len(ptn_seq)
                and ptn_seq[match_start_ind] not in 'NZ'
                and ptn_seq[match_start_ind] == ptn_seq[match_start_ind + 1]):
                return True

        elif slp_type == 'I':
            # Implicit that ptn_seq[match_start_ind + 2] not in 'NZ'
            if (match_start_ind + 2 < len(ptn_seq)
                and ptn_seq[match_start_ind] != 'N'
                and ptn_seq[match_start_ind + 1] not in 'NZ'
                and ptn_seq[match_start_ind] != ptn_seq[match_start_ind + 1]
                and ptn_seq[match_start_ind + 1] == ptn_seq[match_start_ind + 2]):
                return True

        elif slp_type == 'D':
            if (match_start_ind + 1 < len(ptn_seq)
                and 'Z' not in ptn_seq
                and ptn_seq[match_start_ind] not in 'N'
                and ptn_seq[match_start_ind + 1] not in 'N'
                and ptn_seq[match_start_ind] != ptn_seq[match_start_ind + 1]):
                return True

        else:
            raise ValueError('Unimplemented slippage type: {}. {}'.format(slp_type, slp_ptn_tup))

        return False

    

SlippagePtnTuple = collections.namedtuple('SlippagePtnTuple', 
                                          ('ptn_seq', 
                                           'slippage_type', 
                                           'match_start_ind'))

def iterate_slippage_ptn_file(slp_ptn_fn):
    with open(slp_ptn_fn, 'r') as slp_ptn_file:
        for line in slp_ptn_file:
            fields = line.strip().split()
            if len(fields) != 3:
                raise ValueError('Slippage pattern record should have 3 fields. ' + line)

            slp_ptn_tup = SlippagePtnTuple(ptn_seq = fields[0],
                                           slippage_type = fields[1],
                                           match_start_ind = int(fields[2]))

            if not PtnSlippageTbl.is_valid_slp_ptn_tup(slp_ptn_tup):
                raise ValueError('Invalid slippage pattern record: {}'.format(slp_ptn_tup))

            yield slp_ptn_tup

    
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

    arg_parser.add_argument('ptn_tup_iterator', metavar = '<slippage pattern table>', 
                            type = iterate_slippage_ptn_file)

    arg_parser.add_argument('rna_parsed_tbl_iterator', 
                            metavar = '<rna parsed file>', 
                            type = rnautil.iterate_rna_parsed_record,
                            help = 'All DNA/RNA promoter region sequence count file. '
                                   'Tab delimited. '
                                   'Columns from left to right are RNA, DNA, XX, raw count, tn count, and XX. '
                                   'XX denotes deprecated fields. ')

    arg_parser.add_argument('ofile', metavar = '<output file>', type = argparse.FileType('w'))

    args = arg_parser.parse_args()

    ptn_tup_tup = tuple(args.ptn_tup_iterator)
    lib_slp_tbl = LibrarySlippageTable(args.rna_parsed_tbl_iterator, args.beg_seq)
    ptn_slp_tbl = lib_slp_tbl.get_ptn_slippage_tbl(ptn_tup_tup)
    ptn_slp_tbl_str = ptn_slp_tbl.fmt_slippage_tbl()
    args.ofile.write(ptn_slp_tbl_str)
    args.ofile.close()

if __name__ == '__main__':
    main()
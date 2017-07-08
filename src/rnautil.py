import collections

def check_match(seqa, seqb):
    length = min(len(seqa), len(seqb))
    for i in xrange(-1, -length - 1, -1):
        if seqa[i] != seqb[i]:
            return(False)
    return True

def seq_rmatch(seqa, seqb):
    return check_match(seqa, seqb)


# Format of RNA parsed file:
# 0. RNA promoter seq
# 1. DNA promoter seq (fixed length)
# 2. Start position, obsolete
# 3. Raw count
# 4. Digital tag normalized count
# 5. Matched or not, obsolete
# (DNA, RNA) might be duplicated
def get_dna_prmt_seq_len(i_fn):
    i_file = open(i_fn)
    line = i_file.readline()
    i_file.close()

    fields = line.strip().split()
    assert len(fields) == 6
    dna_prmt_seq = fields[1]
    return len(dna_prmt_seq)

def iterate_rna_parsed_file(rp_fn):
    rp_file = open(rp_fn, 'r')
    for line in rp_file:
        fields = line.strip().split()
        assert len(fields) == 6
        yield tuple(fields)

    rp_file.close()

# Require that U in RNA is converted to T throught the analysis. 
def is_valid_seq(seq):
    if len(seq) == 0:
        return False

    for nt in seq:
        if nt not in 'ACGT':
            return False
    return True

RNAParsedRecord = collections.namedtuple('RNAParsedRecord', 
                                         ('rna_prmt_seq', 
                                          'dna_prmt_seq', 
                                          'raw_cnt',
                                          'tn_cnt'))

def iterate_rna_parsed_record(rp_fn):
    rp_file = open(rp_fn, 'r')

    for line in rp_file:
        fields = line.strip().split()
        assert len(fields) == 6
        
        rna_parsed_record = RNAParsedRecord(rna_prmt_seq = fields[0],
                                            dna_prmt_seq = fields[1],
                                            raw_cnt = int(fields[3]),
                                            tn_cnt = int(fields[4]))

        assert is_valid_seq(rna_parsed_record.rna_prmt_seq)
        assert is_valid_seq(rna_parsed_record.dna_prmt_seq)
        assert rna_parsed_record.raw_cnt >= rna_parsed_record.tn_cnt

        yield rna_parsed_record

    rp_file.close()

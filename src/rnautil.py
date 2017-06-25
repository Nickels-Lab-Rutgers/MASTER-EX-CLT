def check_match(seqa, seqb):
    length = min(len(seqa), len(seqb))
    for i in xrange(-1, -length - 1, -1):
        if seqa[i] != seqb[i]:
            return(False)
    return True


# Format of RNA parsed file:
# 0. RNA promoter seq
# 1. DNA promoter seq (fixed length)
# 2. Start position, obsolete
# 3. Raw count
# 4. Digital tag normalized count
# 5. Matched or not, obsolete
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

#!/usr/bin/env python2.7
import argparse
import fastqutil


# For testing, use global variables.
BEG_SEQ = "ATAGAACTTTAGGCACCCCAGGCTTGACACTTTATGCTTCGGCTCGTATAATGTGTGGAA" 
MID_SEQ = "GATAACAATTTCAACAAT"
END_SEQ = ""

def parse_dna_seq(seq, qscore, qcutoff, plen, blen):
    seq = seq.strip()
    seqlen = len(seq)
    qscore = qscore.strip()
    assert seqlen == len(qscore), "seq length (%s) != qscore length (%s)" % (seq, qscore)

    beg_start = seq.find(BEG_SEQ)
    if (beg_start == -1):
        return "BEG_SEQ_NOT_FOUND"

    # beg_start is the length before BEG_SEQ
    if beg_start + len(BEG_SEQ) + plen + len(MID_SEQ) + blen + len(END_SEQ) > seqlen:
        return "SHORT_SEQ"

    # beg-p-mid-b-end
    # pstart is random promoter region start (0-based inc)
    # bstart is barcode start
    pstart = beg_start + len(BEG_SEQ)
    mid_start = pstart + plen
    bstart = mid_start + len(MID_SEQ)
    end_start = bstart + blen
    end_end = end_start + len(END_SEQ)

    seq_ex_mid = seq[mid_start:bstart]
    if seq_ex_mid != MID_SEQ:
        return "MID_SEQ_NOT_FOUND"

    seq_ex_end = seq[end_start:end_start + len(END_SEQ)]
    if seq_ex_end != END_SEQ:
        return "END_SEQ_NOT_FOUND"
    
    if min(qscore[beg_start:end_end]) < qcutoff or "N" in seq[beg_start:end_end]:
        return "QUAL_FAILED"
    
    seq_ex_p = seq[pstart:mid_start]
    seq_ex_b = seq[bstart:end_start]

    return (seq_ex_p, seq_ex_b)

class DNASeuenceContainer:
    def __init__(self):
        self.dna_seq_dict = {}
        self.num_seq = 0
    def insert_seq(self, prmt_seq, brcd_seq):
        self.num_seq += 1
        if brcd_seq not in self.dna_seq_dict:
            self.dna_seq_dict[brcd_seq] = {prmt_seq : 1}
        elif prmt_seq not in self.dna_seq_dict[brcd_seq]:
            self.dna_seq_dict[brcd_seq][prmt_seq] = 1
        else:
            self.dna_seq_dict[brcd_seq][prmt_seq] += 1

    # Output barcode\tpromoter\tcount
    # Return stats
    # cutoff >= 
    def output_seq(self, output_dna_seq_fn, ratio_cutoff, cnt_cutoff):
        num_brcd = 0
        num_low_count_brcd = 0
        num_low_count_reads = 0
        num_multi_map_brcd = 0
        num_multi_map_reads = 0
        num_valid_brcd = 0
        num_valid_reads = 0

        or_file = open(output_dna_seq_fn, 'w')
        prmt_set = set()
        
        for brcd_seq in sorted(self.dna_seq_dict.keys()):
            num_brcd += 1
            prmt_dict = self.dna_seq_dict[brcd_seq]
            prmt_count_list = prmt_dict.values()
            sum_prmt_count = sum(prmt_count_list)
            max_prmt_count = max(prmt_count_list)
            prmt_max_pct = float(max_prmt_count) / sum_prmt_count

            if prmt_max_pct < ratio_cutoff:
                num_multi_map_brcd += 1
                num_multi_map_reads += sum_prmt_count
            elif max_prmt_count < cnt_cutoff:
                num_low_count_brcd += 1
                num_low_count_reads += sum_prmt_count
            else:
                num_valid_brcd += 1
                num_valid_reads += max_prmt_count
                for prmt_seq in prmt_dict:
                    if prmt_dict[prmt_seq] == max_prmt_count:
                        or_file.write('%s\t%s\t%d\n' % (brcd_seq, prmt_seq, max_prmt_count))
                        prmt_set.add(prmt_seq)
                    else:
                        # ISSUE: Diversity should only be incremented once
                        num_multi_map_brcd += 1
                        num_multi_map_reads += prmt_dict[prmt_seq]

        or_file.close()

        stats = "Number of total barcodes: %d\n" % num_brcd
        # Avoid divide by 0
        if num_brcd == 0:
            num_brcd = 1
        stats += "Number of valid barcodes: %d (%s%%)\n" % (
            num_valid_brcd, 
            "{:.4f}".format(float(num_valid_brcd) / num_brcd * 100))
        parsed_read_cnt = self.num_seq
        if parsed_read_cnt == 0:
            parsed_read_cnt = 1
        stats += "Number of valid reads: %d (%s%%)\n" % (
            num_valid_reads, 
            "{:.4f}".format(float(num_valid_reads) / parsed_read_cnt * 100))
        stats += "Number of promoter regions: %d\n" % len(prmt_set)
        stats += "Number of multi-mapped barcodes: %d (%s%%)\n" % (
            num_multi_map_brcd, 
            "{:.4f}".format(float(num_multi_map_brcd) / num_brcd * 100))
        stats += "Number of multi-mapped reads: %d (%s%%)\n" % (
            num_multi_map_reads, 
            "{:.4f}".format(float(num_multi_map_reads) / parsed_read_cnt * 100))
        stats += "Number of low count barcodes: %d (%s%%)\n" % (
            num_low_count_brcd,
            "{:.4f}".format(float(num_low_count_brcd) / num_brcd * 100))
        stats += "Number of low count reads: %d (%s%%)\n" % (
            num_low_count_reads,
            "{:.4f}".format(float(num_low_count_reads) / parsed_read_cnt * 100))
        return stats


def parse_dna_fastq_files(qcutoff, plen, blen, or_fn, os_fn, ifn_list, rcutoff = 10):
    num_total_reads = 0
    num_struct_failed_reads = 0
    num_qual_failed_reads = 0
    num_parsed_reads = 0

    # {barcode : {promoter : count, ... }, ...}
    dna_seq_container = DNASeuenceContainer()
    
    for seq_fn in ifn_list:
        for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(seq_fn):
            num_total_reads += 1

            seq_parse_result = parse_dna_seq(seq, qscore, qcutoff, plen, blen)
            if seq_parse_result in ("BEG_SEQ_NOT_FOUND", "SHORT_SEQ", 
                "MID_SEQ_NOT_FOUND", "END_SEQ_NOT_FOUND"):
                num_struct_failed_reads += 1
            elif seq_parse_result == "QUAL_FAILED":
                num_qual_failed_reads += 1
            else:
                num_parsed_reads += 1
                prmt_seq, brcd_seq = seq_parse_result
                dna_seq_container.insert_seq(prmt_seq, brcd_seq)

    stats = "Analyzed files: " + str(ifn_list) + "\n"
    stats += "Beginning fixed sequence: %s\n" % BEG_SEQ
    stats += "Random promoter region length: %d\n" % plen
    stats += "Middle fixed sequence: %s\n" % MID_SEQ
    stats += "Random barcode region length: %d\n" % blen
    stats += "End fixed sequence: %s\n\n" % END_SEQ
    stats += "Total number of reads: %d\n" % num_total_reads
    # Avoid divide by 0
    if num_total_reads == 0:
        num_total_reads = 1
    stats += "Number of parsed reads: %d (%s%%)\n" % (
        num_parsed_reads, 
        "{:.4f}".format(float(num_parsed_reads) / num_total_reads * 100))
    stats += "Number of structure failed reads: %d (%s%%)\n" % (
        num_struct_failed_reads, 
        "{:.4f}".format(float(num_struct_failed_reads) / num_total_reads * 100))
    stats += "Number of quality failed reads: %d (%s%%)\n\n" % (
        num_qual_failed_reads, 
        "{:.4f}".format(float(num_qual_failed_reads) / num_total_reads * 100))

    stats += "For all %d parsed reads:\n" % num_parsed_reads
    stats += dna_seq_container.output_seq(or_fn, 0.9, rcutoff)

    with open(os_fn, 'w') as os_file:
        os_file.write(stats)

    return num_total_reads



def main():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-q', '--qcutoff', type = int, default = 0,
                            help = 'FASTQ per base Phred quality score cutoff.'
                                   'Discard reads with quality score < qcutoff')

    arg_parser.add_argument('-r', '--rcutoff', type = int, default = 10,
                            help = 'Read count cutoff of barcodes. '
                                   'Discard barcode-promoter pairs with < cutoff '
                                   'reads.')
    
    arg_parser.add_argument('plen', metavar = '<random promoter region length>',
                            type = int)

    arg_parser.add_argument('blen', metavar = '<barcode region length>',
                            type = int)

    arg_parser.add_argument('or_fn', metavar = '<output DNA parsed file name>')

    arg_parser.add_argument('os_fn', metavar = '<output DNA stats file name>')

    arg_parser.add_argument('ifn_list', nargs = '+',
                            metavar = '<FASTQ files (space separated)>')

    args = arg_parser.parse_args()

    parse_dna_fastq_files(fastqutil.phread_quality_int_to_char(args.qcutoff),
                          args.plen, args.blen, args.or_fn, args.os_fn,
                          args.ifn_list, rcutoff = args.rcutoff)

    return 0

if __name__ == '__main__':
    main()

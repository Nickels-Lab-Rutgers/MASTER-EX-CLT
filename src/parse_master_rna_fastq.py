#!/usr/bin/env python2.7
import argparse
import fastqutil
import rnautil


# For testing, use global variables.
MID_SEQ = "TTCAACAAT"
END_SEQ = "TGGAATTCTCGGGTG"

def parse_rna_seq(seq, qscore, qcutoff, dlen, plen, blen):
    seq = seq.strip()
    seqlen = len(seq)
    qscore = qscore.strip()
    assert seqlen == len(qscore), "seq length (%s) != qscore length (%s)" % (seq, qscore)

    # Assume that digital tag always exist, so find MID_SEQ from dlen + 1.
    # dlen + 1 because at least 1 base of promoter region needs to be identified.

    # Starting index to search for mid_seq. Inclusive. 
    # So it is also the length to skip at the beginning.
    mid_search_start_ind = dlen + 1

    # Traverse all occurence of mid_seq in the read. Return the first valid one.
    # Assume that the end_seq is long enough to exclude false positive ones. 
    while True:
        mid_start = seq.find(MID_SEQ, mid_search_start_ind)
        if (mid_start == -1):
            return "MID_SEQ_NOT_FOUND"
    
        # digital_tag-p-mid-b-end
        if mid_start + len(MID_SEQ) + blen + len(END_SEQ) > seqlen:
            return "SHORT_SEQ"
        
        # pstart is random promoter region start (0-based inc)
        # bstart is barcode start
        pstart = dlen
        bstart = mid_start + len(MID_SEQ)
        end_start = bstart + blen
        end_end = end_start + len(END_SEQ)
    
        seq_ex_end = seq[end_start:end_end]
        if seq_ex_end != END_SEQ:
            #return "END_SEQ_NOT_FOUND"
            mid_search_start_ind = mid_start + 1
        else:
            break
    
    if min(qscore[0:end_end]) < qcutoff or "N" in seq[0:end_end]:
        return "QUAL_FAILED"
    
    seq_ex_dt = seq[0:dlen]
    seq_ex_p = seq[pstart:mid_start]
    seq_ex_b = seq[bstart:end_start]

    return (seq_ex_dt, seq_ex_p, seq_ex_b)

class RNASequenceContainer:
    def __init__(self, dna_parsed_tpl_fn):
        self.rna_seq_dict = {}
        self.dna_tpl_dict = self.get_dna_bpdict(dna_parsed_tpl_fn)
        self.num_seq = 0
    
    @staticmethod
    def get_dna_bpdict(idna_parsed_fn):
        # {barcode : (promoter, count), ...}
        dna_bpdict = {}
        with open(idna_parsed_fn, 'r') as ifile:
            for line in ifile:
                fields = line.strip().split()
                assert len(fields) == 3
                brcd_seq = fields[0]
                prmt_seq = fields[1]
                prmt_count = int(fields[2])
                assert brcd_seq not in dna_bpdict
                dna_bpdict[brcd_seq] = (prmt_seq, prmt_count)
        return dna_bpdict

    def insert_seq(self, dtag_seq, prmt_seq, brcd_seq):
        self.num_seq += 1
        if brcd_seq not in self.rna_seq_dict:
            self.rna_seq_dict[brcd_seq] = {prmt_seq : [1, {dtag_seq : 1}]}
        elif prmt_seq not in self.rna_seq_dict[brcd_seq]:
            self.rna_seq_dict[brcd_seq][prmt_seq] = [1, {dtag_seq : 1}]
        elif dtag_seq not in self.rna_seq_dict[brcd_seq][prmt_seq][1]:
            self.rna_seq_dict[brcd_seq][prmt_seq][0] += 1
            self.rna_seq_dict[brcd_seq][prmt_seq][1][dtag_seq] = 1
        else:
            self.rna_seq_dict[brcd_seq][prmt_seq][0] += 1
            self.rna_seq_dict[brcd_seq][prmt_seq][1][dtag_seq] += 1

    # Write saved RNA sequences
    # Return RNA stats
    def output_seq(self, output_fn):
        num_brcd = 0
        num_brcd_notpl = 0
        num_reads_notpl = 0
        sum_dtag_count_notpl = 0
        num_valid_brcd = 0
        num_valid_reads = 0
        sum_dtag_count_valid = 0
    
        identified_dna_prmt_set = set()

        or_file = open(output_fn, 'w')
        for brcd_seq, prmt_dict in self.rna_seq_dict.iteritems():
            num_brcd += 1

            if brcd_seq not in self.dna_tpl_dict:
                num_brcd_notpl += 1
                for prmt_seq in prmt_dict:
                    prmt_count_list = prmt_dict[prmt_seq]
                    num_reads_notpl += prmt_count_list[0]
                    sum_dtag_count_notpl += len(prmt_count_list[1])
            else:
                dna_prmt_seq = self.dna_tpl_dict[brcd_seq][0]
                identified_dna_prmt_set.add(dna_prmt_seq)
                num_valid_brcd += 1
                for prmt_seq, prmt_count_list in prmt_dict.iteritems():
                    if rnautil.seq_rmatch(prmt_seq, dna_prmt_seq):
                        prmt_rmatch = 1
                    else:
                        prmt_rmatch = 0
                    start_ind = len(dna_prmt_seq) - len(prmt_seq) + 1

                    prmt_read_count = prmt_count_list[0]
                    prmt_dtag_conut = len(prmt_count_list[1])

                    num_valid_reads += prmt_read_count
                    sum_dtag_count_valid += prmt_dtag_conut
                    or_file.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (prmt_seq, 
                        dna_prmt_seq, start_ind, prmt_read_count, prmt_dtag_conut,
                        prmt_rmatch))
        or_file.close()
        
        sum_dtag_counts = sum_dtag_count_valid + sum_dtag_count_notpl

        stats = "Sum of digital tag counts: %d\n" % sum_dtag_counts
        # Avoid divide by 0
        if sum_dtag_counts == 0:
            sum_dtag_counts = 1

        stats += "Number of total barcodes: %d\n" % num_brcd
        if num_brcd == 0:
            num_brcd = 1

        stats += "Number of DNA template promoter regions: %d\n" % len(identified_dna_prmt_set)
        stats += "Number of valid barcodes: %d (%s%%)\n" % (
            num_valid_brcd, 
            "{:.4f}".format(float(num_valid_brcd) / num_brcd * 100))
        seq_cnt = self.num_seq
        if seq_cnt == 0:
            seq_cnt = 1
        stats += "Number of valid reads: %d (%s%%)\n" % (
            num_valid_reads, 
            "{:.4f}".format(float(num_valid_reads) / seq_cnt * 100))
        stats += "Sum of valid digital tag counts: %d (%s%%)\n" % (
            sum_dtag_count_valid, 
            "{:.4f}".format(float(sum_dtag_count_valid) / sum_dtag_counts * 100))

        stats += "Number of no-dna-template barcodes: %d (%s%%)\n" % (
            num_brcd_notpl, 
            "{:.4f}".format(float(num_brcd_notpl) / num_brcd * 100))
        stats += "Number of no-dna-template reads: %d (%s%%)\n" % (
            num_reads_notpl, 
            "{:.4f}".format(float(num_reads_notpl) / seq_cnt * 100))
        stats += "Sum of no-dna-template digital tag counts: %d (%s%%)\n" % (
            sum_dtag_count_notpl, 
            "{:.4f}".format(float(sum_dtag_count_notpl) / sum_dtag_counts * 100))
        return stats



def parse_rna_fastq_files(qcutoff, dlen, plen, blen, idna_parsed_fn, or_fn, os_fn, ifn_list):
    num_total_reads = 0
    num_struct_failed_reads = 0
    num_qual_failed_reads = 0
    num_parsed_reads = 0

    # {brcd : {prmt : [raw_cnt, {dtag : rawcnt, ...}], ...}, ...}
    rna_seq_container = RNASequenceContainer(idna_parsed_fn)

    for seq_fn in ifn_list:
        for seqid, seq, rseqid, qscore in fastqutil.iterate_fastq_file(seq_fn):
            num_total_reads += 1

            seq_parse_result = parse_rna_seq(seq, qscore, qcutoff, dlen, plen, blen)
            if seq_parse_result in ("SHORT_SEQ", "MID_SEQ_NOT_FOUND", 
                "END_SEQ_NOT_FOUND"):
                num_struct_failed_reads += 1
            elif seq_parse_result == "QUAL_FAILED":
                num_qual_failed_reads += 1
            else:
                num_parsed_reads += 1
                dtag_seq, prmt_seq, brcd_seq = seq_parse_result
                rna_seq_container.insert_seq(dtag_seq, prmt_seq, brcd_seq)

    # Generate read stats 
    stats = "RNA files: " + str(ifn_list) + "\n"
    stats += "DNA template file: " + idna_parsed_fn + "\n"
    stats += "Digital tag length: %d\n" % dlen
    stats += "Promoter region length: %d\n" % plen
    stats += "Middle sequence: " + MID_SEQ + "\n"
    stats += "Barcode region length: %d\n" % blen
    stats += "End sequence: " + END_SEQ + "\n"
    stats += "Quality score cutoff: %d\n\n" % (ord(qcutoff) - 33)
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
    stats += rna_seq_container.output_seq(or_fn)

    with open(os_fn, 'w') as os_file:
        os_file.write(stats)

    return num_total_reads

def main():
    arg_parser = argparse.ArgumentParser()

    arg_parser.add_argument('-q', '--qcutoff', type = int, default = 0,
                            help = 'FASTQ per base Phred quality score cutoff.'
                                   'Discard reads with quality score < qcutoff')

    arg_parser.add_argument('dlen', metavar = '<digital tag length>',
                            type = int)

    arg_parser.add_argument('plen', metavar = '<random promoter region length>',
                            type = int)

    arg_parser.add_argument('blen', metavar = '<barcode region length>',
                            type = int)

    arg_parser.add_argument('idna_parsed_fn', metavar = '<input DNA parsed file>')

    arg_parser.add_argument('or_fn', metavar = '<output RNA parsed file name>')

    arg_parser.add_argument('os_fn', metavar = '<output RNA stats file name>')

    arg_parser.add_argument('ifn_list', nargs = '+',
                            metavar = '<FASTQ files (space separated)>')

    args = arg_parser.parse_args()
    
    parse_rna_fastq_files(fastqutil.phread_quality_int_to_char(args.qcutoff),
                          args.dlen, args.plen, args.blen, args.idna_parsed_fn,
                          args.or_fn, args.os_fn, args.ifn_list)

    return 0

if __name__ == '__main__':
    main()

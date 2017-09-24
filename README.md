# MASTER Extended Command Line Tools
This package includes a set of command line tools for analyzing MAssively Systematic Transcript End Readout (MASTER) data. The purpose is to determine the number of RNA reads emanating from each position of DNA template. 

Two main steps of general MASTER data analysis are DNA template analysis and 5' RNA-Seq analysis. 

## DNA template analysis
The purpose of analyzing the sequencing results of DNA template library is to associate the 7-bp randomized TSS-region with a corresponding randomized 15-bp barcode region. 

`parse_master_dna_fastq.py` is used to analyze sequencing results of DNA template libraries. The output TSS-regions and barcodes file is tab delimited file with barcode, TSS-region, and count written from left to right. This program will also output a stats file of analyzed DNA template library. 

`count_parsed_dna.py` is used to count the number of reads of each TSS-region in DNA parsed results. 

## 5' RNA-Seq analysis
The purpose of analyzing 5' RNA-Seq results is to identify the DNA templates and transcription start position of RNA reads.

`parse_master_rna_fastq.py` is used to analyze 5' RNA-Seq sequencing results. 

`ex_count_parsed_rna.py` is used to count the number of RNA reads emanating from each position of DNA template. 

## Other command line tools

`attach_dna_prmt_suffix.py`
Attach suffix sequence to the TSS-region sequence in DNA parsed results. 

`count_seq_prefix.py`
Count the number of reads with different prefix sequence in FASTQ file.

`count_spike_in_seq.py`
Count the number of spike-in sequence reads in FASTQ file. 

`gen_merged_rna_prmt_seq_tbl.py`
Generate the 5' RNA-Seq read summary files using parsed results. 

`prefix_grep_fastq.py`
Select sequencing reads in FASTQ file with specific prefix. 

`profile_slippage.py`
Profile the transcription slippage in 5' RNA-Seq reads. 
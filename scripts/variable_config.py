#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generated Friday, September 24, 2021: 10:49:24
"""

tntblast_location = ""
assay_type = 'PCR' # str, Assay format to simulate
min_primer_Tm = 40 # int, Minimum primer Tm
min_probe_Tm = 40 # int, Minimum probe Tm
max_amplicon_length = 1000 # int, Maximum amplicon length
minimum_target_length = 29000 # int, Minimum length of target sequences. Doesn't filter GISAID
best_match = True # Only save the best match, in Tm, between a query and target



# Paths
resource_directory =  "/home/adanm/indiv_eval/parse_test/resources"
results_directory =  "/home/adanm/indiv_eval/parse_test/Results"
tnt_results_directory =  "/home/adanm/indiv_eval/parse_test/TNT_Results"
fasta_directory =  "/panfs/biopan04/scratch-adanm/sars_cov_2/test_gisaid"

# Files
assay_file =  "/home/adanm/indiv_eval/parse_test/resources/assays.txt"
sequence_file =  "/panfs/biopan04/scratch-adanm/sars_cov_2/test_gisaid/sequences.fasta"

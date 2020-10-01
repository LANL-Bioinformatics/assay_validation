#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:12:38 2020

This script runs the assay monitor functionality for this package.  Using the
    fasta files downloaded by am_download.py and certain resource files (detailed
    below), assays (provided in a resource file) will be evaluated against each
    genome (previously downloaded fastas) using ThermonucleotideBLAST.  Results
    will be output as full results in Assay_Results.json, as a summary in
    summary_table.json, and as cross referenced results in match_table.csv.
    Summary stats are also output to db_stats.json and db_totals.json.  These
    files are used downstream in this package for producing visualizations.


This script calls ThermonucleotideBLAST (TNTBLAST) which was created by Jason
    Gans at Los Alamos National Laboratory with the BSD 3-Clause License.
    Please contact Jason at jgans@lanl.gov

    The reference for ThermonucleotideBLAST is: "Improved assay-dependent
    searching of nucleic acid sequence databases." by J. D. Gans and M.
    Wolinsky, Nucleic Acids Res. 2008 Jul;36(12):e74. doi: 10.1093/nar/gkn301.

    Much thanks to Jason Gans for creating the algorithm that is the core
    functionality for this package.

    To download and install TNTBLAST, please go to
    https://github.com/jgans/thermonucleotideBLAST


Dependencies:
    db_stats.py
    new_tnt_parse.py
    assay_sum_table.py
    seq_db_filter_funx.py
    TNTBLAST


The script requires the following arguments:

        '-r', '--resource_directory' <resource_directory>
            This is the directory which must contain certain necessary files
            and to which other resource files will be created.
            Use the format:  /XXX/XXX/XXX/

        '-R', '--results_directory' <results_directory>
            This is the directory to which results files will be created.
            Use the format:  /XXX/XXX/XXX/


        '-f', '--fasta_directory' <fasta_directory>
            This is the directory which contains the fna sequence files.
            Use the format:  /XXX/XXX/XXX/

        '-t', '--tnt_results_directory' <tnt_results_directory>
            This is the directory to which results files will be created when
            TNTBLAST is run.
            Use the format:  /XXX/XXX/XXX/



These files must be present in the resource directory:
    "assays.txt"
        A list of assays with corresponding oligo sequences.  Each to a line.
        Order of Oligo sequences: Forward primer Reverse primer, Probe.  Ex:
CDC-2019-nCoV_N1 GACCCCAAAATCAGCGAAAT TCTGGTTACTGCCAGTTGAATCTG ACCCCGCATTACGTTTGGTGGACC

    "reduced_assays.txt"
        A reduced list of assays to be used in the final results.  Same format
        as "assays.txt"

    "del_ct_table.txt"
        A table of delta Ct values from { Li B, Kadura I, Fu DJ, Watson DE.
        Genotyping with TaqMAMA. Genomics. 2004 Feb 1;83(2):311-20. }.
        Tab separated.  Headers for Rows and columns.  First entry: "Row"
        Example first two lines:
Row     CC      GC      AC      TC      CG      GG      AG      TG      CA      GA      AA      TA      CT      GT      AT      TT
CC      0.0     0.3     0.5     -0.3    6.3     17.5    19.2    11.9    5.1     15.7    11.4    12.4    0.4     11.0    10.4    3.7

@author: adanm
"""
import re
import os
import sys
import json
import time
import db_stats
import datetime
import argparse as ap
from sys import platform
import new_tnt_parse as nuparse
import assay_sum_table as assum
import seq_db_filter_funx as filt_db



''' 'Global' Variables '''
# Filter Options
minimum_file_length_considered = 29000

# Dictionary Labels
assay = "Assay"
assay_name = "Name"
FP = "Forward Primer"
RP = "Reverse Primer"
Pr = "Probe"
gp = "gaps"
mm = "mismatches"

# TNTBLAST Options
assay_type = "PCR"
temperature = 40
amp_length = 1000



''' Methods '''

def generate_fna_lists(fasta_dir):

    ''' Get list of all genbank files '''
    genbank_file_list = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir)
                 if re.match(r'.*\.fna', f)
                 and not re.match(r'EPI.*\.fna', f)]
    genbank_file_list.sort()

    return genbank_file_list


def run_tnt_blast(sequence_file, results_directory, tntblast_location, assay_file):
    seq_file_prefix = sequence_file.split(".")[0]
    seq_acc = seq_file_prefix.split("/")[-1]
    tnt_results_file = "{}_results.out".format(seq_acc)
    tnt_results_path = os.path.join(results_directory, tnt_results_file)
    tnt_runline = "{}tntblast -A {} -i {} -d {} -e {} -E {} -l {} -o {} --best-match >> std.err".format(
            tntblast_location, assay_type, assay_file, sequence_file, temperature,
            temperature, amp_length, tnt_results_path)
    os.system("{}".format(tnt_runline))

    return tnt_results_path




def get_arguments():
    parser = ap.ArgumentParser(prog='assay_monitor.py', description="""Assay Monitor""")
    parser.add_argument('-r', '--resource_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="resource files directory")
    parser.add_argument('-R', '--results_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="results files directory")
    parser.add_argument('-f', '--fasta_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="location of fasta files")
    parser.add_argument('-t', '--tnt_results_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="location of tnt results")
    return parser.parse_args()



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    big_tic = time.perf_counter()

    ''' Get command line arguments '''
    args = get_arguments()

    resource_dir = args.resource_directory[0]
    results_dir = args.results_directory[0]
    fasta_dir = args.fasta_directory[0]
    tnt_res_dir = args.tnt_results_directory[0]

    # TNTBLAST location
    if platform == "darwin":
        tntblast_location = "/Users/adanm/bin/"     # On Mac
    elif platform == "linux":
        tntblast_location = ""                      # On Rustang

    # File Names
    db_stats_file = os.path.join(resource_dir, "db_stats.json")
    db_totals_file = os.path.join(resource_dir, "db_totals.json")
    seq_file = os.path.join(resource_dir, "seq_file_list.txt")
    meta_filename = os.path.join(resource_dir, "nextstrain_ncov_metadata.txt")
    sum_data_file = os.path.join(resource_dir, "esummary_ncbi.txt")
    assay_file = os.path.join(resource_dir, "assays.txt")
    short_assay_list_file = os.path.join(resource_dir, "reduced_assays.txt")
    output_file = os.path.join(results_dir, "summary_table.json")
    results_file = os.path.join(results_dir, "Assay_Results.json")
    combined_fasta = os.path.join(fasta_dir, "All_Seqs.fasta")
    three_prime_table = os.path.join(resource_dir, "del_ct_table.txt")
    negatives_list_output_file = os.path.join(resource_dir, "assay_failures.tsv")

    # Variable
    del_ct_threshhold = 2.0


    tic = time.perf_counter()
    ''' ####################### Generate Databases ######################### '''

    genbank_file_list = generate_fna_lists(fasta_dir)

    ''' ######################## Filter Databases ########################## '''

    file_list, genbank_tallies, total_list_len = filt_db.filter_dbs_output_all(
        genbank_file_list, minimum_file_length_considered)

    ''' ##################### Combine Sequence Files ####################### '''

    ''' Combine fna sequence files into single fasta file '''
    print("\nCombining target fna files from database into one fasta...")
    accession_list = []
    with open(combined_fasta, 'w') as write_combine:
        for fna_file in file_list:
            accession = fna_file.split('/')[-1].split('.')[0]

            with open(fna_file, 'r') as read_fna:
                for line in read_fna.readlines():
                    print(line.strip(), file=write_combine)
            accession_list.append(accession)

    ''' ##################### Output DB Stats Files ######################## '''
    db_stats.generate_db_stats(genbank_tallies, total_list_len, db_stats_file, db_totals_file)

    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Generate DB, Filter DB, Combining Files, and DB Stats Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))






    ''' ####################### TNT_BLAST Testing ######################### '''
    tic = time.perf_counter()

    ''' Initialize multilevel dictionary for tnt assay results '''
    Full_Dict = {assay : []}


    ''' Run TNT_BLAST '''
    print("\nRunning tntblast against all sequences in combined file...")
    tnt_results_list = []
    tnt_results = run_tnt_blast(combined_fasta, tnt_res_dir, tntblast_location, assay_file)
    tnt_results_list.append(tnt_results)

    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("TNT_BLAST Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))





    ''' Parse tnt output for True Positives into multi-level dictionary for json file: '''
    tic = time.perf_counter()

    print("\nParsing TNT Output into Metadata Dictionary...")
    for tnt_result in tnt_results_list:
        Full_Dict = nuparse.parse_tnt_results(Full_Dict, tnt_result)


    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Parse Output Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))





    ''' False Negatives '''
    ''' Parse and Print false negative results '''
    tic = time.perf_counter()

    print("\nDetermining Negatives and Adding To Metadata Dictionary...")
    for accession in accession_list:
        Full_Dict = nuparse.determine_negatives(assay_file, accession, Full_Dict, "True Positives", "False Negatives")


    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Determine Negatives Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))




    ''' Three Prime Filter '''
    ''' Move 3' MMs (above del_Ct threshold) from TP to FN '''
    tic = time.perf_counter()

    print("\nRemoving 3' mismatches over threshold from True Positives...")
    Full_Dict = nuparse.filter_three_prime(Full_Dict, three_prime_table, del_ct_threshhold)


    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Three Prime Filter Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))




    ''' ######################### Output Results ########################## '''

    ''' Get current time for time stamp and filename '''
    current_time = datetime.datetime.now()

    ''' Add timestamp '''
    current_time_string = str(current_time)
    Full_Dict.update({ "Timestamp" : current_time_string })

    ''' Output results as json file '''
    print("\nWriting results to", results_file)
    with open(results_file, 'w') as file_handle:
        print(json.dumps(Full_Dict, sort_keys=True, indent=4), file=file_handle)
    print("Wrote file:", results_file)





    ''' Prep sequence list '''
    nuparse.make_sequence_list_file(seq_file, fasta_dir)


    ''' Create Match Table (Reduced assay list) '''
    tic = time.perf_counter()

    print("\nCreating match table from source:", results_file)
    nuparse.create_match_tables_reduced(
            results_file, results_dir, seq_file, short_assay_list_file)

    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Match Table Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))





    ''' Initialize table dictionary '''
    table_dict = { "data": [ ] }

    ''' Get short assay list '''
    short_assay_list = []
    with open(short_assay_list_file, 'r') as read_file:
        for line in read_file.readlines():
            fields = line.split()
            assay_label = fields[0]
            short_assay_list.append(assay_label)


    ''' Create Summary table '''
    tic = time.perf_counter()

    print("\nCreating summary table json from source:", results_file)
    table_dict = assum.generate_table(results_file, meta_filename, sum_data_file, table_dict, short_assay_list)

    ''' Output results as json file '''
    with open(output_file, 'w') as file_handle:
        print(json.dumps(table_dict, sort_keys=True, indent=4), file=file_handle)
    print("Wrote file:", output_file)

    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Summary Table Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))



    ''' Create Negatives List '''
    tic = time.perf_counter()

    print("\nCreating Negatives List from source:", results_file)
    nuparse.make_negatives_list(results_file, negatives_list_output_file, fasta_dir)

    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("Negatives List Runtime: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))





    big_toc = time.perf_counter()
    elapsed_time = big_toc - big_tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("\nFinished")
    print("Total Elapsed Time: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))


if __name__ == "__main__":
    main()

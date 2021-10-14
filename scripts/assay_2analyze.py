#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:12:38 2020

This script runs the assay monitor functionality for this package.  
    Using the fasta files downloaded by am_download.py and certain 
    resource files (detailed below), assays (provided in a resource 
    file) will be evaluated against each genome (previously downloaded 
    fastas) using ThermonucleotideBLAST.  Results will be output as full
    results in Assay_Results.json, as a summary in summary_table.json, 
    and as cross referenced results in match_table.csv.

This script calls ThermonucleotideBLAST (TNTBLAST) which was created by 
    Jason Gans at Los Alamos National Laboratory with the BSD 3-Clause 
    License.  Please contact Jason at jgans@lanl.gov

    The reference for ThermonucleotideBLAST is: 
    "Improved assay-dependent searching of nucleic acid sequence 
    databases." by J. D. Gans and M. Wolinsky, Nucleic Acids Res. 2008 
    Jul;36(12):e74. doi: 10.1093/nar/gkn301.

    Much thanks to Jason Gans for creating the algorithm that is the 
    core functionality for this package.

    To download and install TNTBLAST, please go to
    https://github.com/jgans/thermonucleotideBLAST


Dependencies:
    aux_funx
    assay_sum_table
    newer_tnt_parse_oldtnt
    TNTBLAST


The script requires the following arguments:

        '-r', '--resource_directory' <resource_directory>
            This is the directory which must contain certain necessary 
            files and to which other resource files will be created.
            Use the format:  /XXX/XXX/XXX/

        '-R', '--results_directory' <results_directory>
            This is the directory where results files will be created.
            Use the format:  /XXX/XXX/XXX/


        '-f', '--fasta_directory' <fasta_directory>
            This is the directory which contains the fna sequence files.
            Use the format:  /XXX/XXX/XXX/

        '-t', '--tnt_results_directory' <tnt_results_directory>
            This is the directory to which results files will be created 
            when TNTBLAST is run.
            Use the format:  /XXX/XXX/XXX/



These files must be present in the resource directory:
    "assays.txt"
        A list of assays with corresponding oligo sequences.  Each to a 
        line.  Order of Oligo sequences: Forward primer Reverse primer, 
        Probe.  Ex:
CDC-2019-nCoV_N1 GACCCCAAAATCAGCGAAAT TCTGGTTACTGCCAGTTGAATCTG \
    ACCCCGCATTACGTTTGGTGGACC

    "reduced_assays.txt"
        A reduced list of assays to be used in the final results.  Same 
        format as "assays.txt"

    "del_ct_table.txt"
        A table of delta Ct values from the reference:
        { Li B, Kadura I, Fu DJ, Watson DE. Genotyping with TaqMAMA. 
        Genomics. 2004 Feb 1;83(2):311-20. }.
        This file should be tab separated.  Headers for Rows and 
        columns.  First entry: "Row"
        Example first two lines:
Row     CC      GC      AC      TC      CG      GG      AG      \
    TG      CA     GA      AA      TA      CT      GT      AT      TT
CC      0.0     0.3     0.5     -0.3    6.3     17.5    19.2    \
    11.9    5.1    15.7    11.4    12.4    0.4     11.0    10.4    3.7

@author: adanm
"""
import os
import sys
import time
import aux_funx
import ujson as json
import argparse as ap
import assay_sum_table as assum
import newer_tnt_parse_oldtnt as nuparse


''' Methods '''

def write_dict_to_json(Full_Dict, results_file):
    print("\nWriting results to", results_file)
    with open(results_file, 'w') as file_handle:
        print(
            json.dumps(
                Full_Dict, sort_keys=True, indent=4
            ), file=file_handle
        )
    print("Wrote file:", results_file)


def print_time(tic, msg):
    toc = time.perf_counter()
    elapsed_time = toc-tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    secs = elapsed_time - minutes * 60 - hours * 3600
    print(
        f"{msg} Runtime: {hours} hours, {minutes} minutes, {secs} seconds"
    )
    print()


def create_filename_namespace(args):
    filenames_ns = ap.Namespace()

    resource_dir = args.resource_directory
    results_dir = args.results_directory
    fasta_directory = args.fasta_directory
    tnt_dir = args.tnt_results_directory
    
    # Filenames
    gisaid_filename = "sequences.fasta"
    genbank_filename = "All_Seqs.fasta"
    if args.use_gisaid:
        sequence_filename = gisaid_filename
        tnt_resultname = "sequences_results.out"
    else:
        sequence_filename = genbank_filename
        tnt_resultname = "All_Seqs_results.out"
    meta_filename = "metadata.tsv"
    accession_filename = "accessions.txt"
    results_filename = "Assay_Results.json"
    three_prime_filename = "del_ct_table.txt"
    sum_table_name = "summary_table.json"
    accession_list_filename = "accession_list.txt"

    # Paths
    filenames_ns.resource_dir = args.resource_directory
    filenames_ns.results_dir = args.results_directory
    filenames_ns.fasta_directory = args.fasta_directory
    filenames_ns.tnt_dir = args.tnt_results_directory
        
    # Files
    filenames_ns.metafile = os.path.join(
            resource_dir, meta_filename)
    filenames_ns.acclist_file = os.path.join(
            results_dir, accession_list_filename)
    filenames_ns.sequence_file = os.path.join(
            fasta_directory, sequence_filename)
    filenames_ns.tnt_result = os.path.join(tnt_dir, tnt_resultname)
    filenames_ns.accession_file = os.path.join(
            results_dir, accession_filename)
    filenames_ns.results_json = os.path.join(
            results_dir, results_filename)
    filenames_ns.three_prime_table = os.path.join(
            resource_dir, three_prime_filename)
    filenames_ns.sum_table_file = os.path.join(
            results_dir, sum_table_name)
    
    # GISAID flag
    filenames_ns.use_gisaid = args.use_gisaid

    return filenames_ns


def get_arguments():
    parser = ap.ArgumentParser(prog=sys.argv[0].split('/')[-1])

    parser.add_argument('-G', '--use_gisaid', action='store_true',
                        help="Use GISAID metadata file.  If not " +
                        "present or using GenBank, omit flag.")
    parser.add_argument('-r', '--resource_directory', metavar='[STR]', 
                        type=str,required=True, 
                        help="resource files directory")
    parser.add_argument('-R', '--results_directory', metavar='[STR]', 
                        type=str, required=True, 
                        help="results files directory")
    parser.add_argument('-f', '--fasta_directory', metavar='[STR]', 
                        type=str, required=True, 
                        help="location of fasta files")
    parser.add_argument('-t', '--tnt_results_directory', 
                        metavar='[STR]', type=str, required=True, 
                        help="location of tnt results")
    parser.add_argument('-p', '--procnum',
                        metavar='[INT]', type=int, required=True,
                        help="Number of processors to use for parallel " +
                        "operations (in three prime filter and summary " +
                        "table generation).")
    parser.add_argument('-d', '--del_ct_threshhold',
                        metavar='[FLOAT]', type=float, default=2.0,
                        help="Delta Ct value to use as threshold in " +
                        "determining thermo mismatches in the three " +
                        "prime filter. Default value is 2.0.")
    return parser.parse_args()



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    big_tic = time.perf_counter()

    # Get command line arguments and generate filename namespace
    args = get_arguments()
    file_ns = create_filename_namespace(args)

    # Variables
    del_ct_threshhold = args.del_ct_threshhold
    num_of_procs = args.procnum
    if args.use_gisaid:
        iso_idx = 0
    else:
        iso_idx = 1


#    ''' ################### TNT_BLAST Parsing ###################### '''
#    tic = time.perf_counter()
#    
#    # Prepare accession list
#    iso_acc_dict = nuparse.get_acc_list(file_ns)
#
#    # Parse tnt output into multi-level dict of true positives
#    print("\nParsing TNT Output into Metadata Dictionary...")
#    Full_Dict = nuparse.brute_parse_tnt_results(
#            file_ns.tnt_result, iso_acc_dict, iso_idx)
#    # Full_Dict = nuparse.parse_tnt_results(
#    #         file_ns.tnt_result, iso_acc_dict, iso_idx)
#
#    # Output results as json file
#    write_dict_to_json(Full_Dict, file_ns.results_json)
#
#    print_time(tic, "Parsing")
#    ''' ################# End TNT_BLAST Parsing #################### '''
#
#
#
#
#    ''' Accession Make '''
#    tic = time.perf_counter()
#    aux_funx.generate_accession_file(file_ns, iso_idx)
#    print_time(tic, "Accession Listmake")
#
#
#    
#    ''' ############### Determine False Negatives ################## '''
#    # Get accession list from accession file
#    accession_list = []
#    with open(file_ns.accession_file) as acc_file:
#        for line in acc_file.readlines():
#            accession_list.append(line.strip())
#    
#    # Get assay monitor data from json file
#    with open(file_ns.results_json) as json_file:
#        Full_Dict = json.load(json_file)
#
#    # Parse and Determine false negative results
#    tic = time.perf_counter()
#    Full_Dict = nuparse.determine_negatives(
#            accession_list, Full_Dict, 
#            "True Positives", "False Negatives")
#    print_time(tic, "Determine Negatives")
#
#    # Output results as json file
#    tic = time.perf_counter()
#    write_dict_to_json(Full_Dict, file_ns.results_json)
#    print_time(tic, "Writing")
#
#    ''' ################## End Determine Negatives ################# '''



    
    ''' #################### Three Prime Filter #################### '''
    tic = time.perf_counter()

    # Get assay monitor data from json file
    with open(file_ns.results_json) as json_file:
        Full_Dict = json.load(json_file)
        
    # Move 3' MMs (above del_Ct threshold) from TP to FN (Parallel)
    Full_Dict = nuparse.set_queues_for_three_prime(
            Full_Dict, file_ns.three_prime_table, 
            del_ct_threshhold, num_of_procs)
    # # Move 3' MMs (above del_Ct threshold) from TP to FN (single thread)
    # Full_Dict = nuparse.filter_three_prime(
    #         Full_Dict, file_ns.three_prime_table, del_ct_threshhold)

    # Output results as json file
    write_dict_to_json(Full_Dict, file_ns.results_json)

    print_time(tic, "Three Prime Filter")
    ''' ################# End Three Prime Filter ################### '''







    ''' ################### Create Summary Table ################### '''
    tic = time.perf_counter()

    # Get assay monitor data from json file
    with open(file_ns.results_json, 'r') as json_file:
        Full_Dict = json.load(json_file)

    # Generate dictionary containing summary results
    # table_dict = assum.generate_table(Full_Dict)
    table_dict = assum.generate_table_parallel(Full_Dict, num_of_procs)

    # Output results as json file
    write_dict_to_json(table_dict, file_ns.sum_table_file)

    print_time(tic, "Summary Table")
    ''' ################# End Create Summary Table ################# '''




    big_toc = time.perf_counter()
    elapsed_time = big_toc - big_tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("\nFinished")
    print(
        "Total Elapsed Time: {} hours, {} minutes, {} seconds".format(
            hours, minutes, seconds))


if __name__ == "__main__":
    main()

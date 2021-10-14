#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:03:28 2021

This script generates files necessary for the pipeline-downstream module 
    (run_tnt.py in the tnt_wrapper repo) to run.
    
    It takes a number of arguments that include full-path locations of the
    necessary files and arguments that specify configuration options for 
    TNTBLAST (which is run by run_tnt.py).  The configuration file is written
    in the indicated directory and will include the variables specified by
    these arguments.  The action with regard to preparation of the sequence
    file depends on the argument chosen.  If "--use_gisaid_data" is indicated
    than it will be assumed that there exists a single fasta file containing
    all the sequences in the location indicated.  If "--use_genbank_data" is
    indicated then it will be assumed that multiple fasta files (with .fna
    affix) for each individual sequence will be present in the indicated 
    directory.  In this case, these files will be assembled into a single
    large fasta.  If the "--use_tax_tree" option is indicated, it will be 
    assumed that there is a directory-based taxonomic tree with its root at
    the location indicated.  The reponse to this has not yet been implemented.
    
    In the case of "--use_gisaid_data", sequences will be filtered to remove
    any that contain the words "bat" or "pangolin" in the header.  Any that
    have a sequence length less than that indicated by "--minimum_target_length"
    will also be filtered out.  The numbers filtered or not will be tallied
    and output as the summary stats files: db_stats.json and db_totals.json.
    These files are used downstream in assay_validation for producing 
    visualizations.
    
Dependencies:
    global_txt.py
    db_stats.py
    seq_db_filter_funx.py
    config_prep.py
    gbk_file_prep.py
    
    
This script requires the following arguments:

    (-c [STR], --config_save_loc [STR])
            directory in which to save written config file
            
    (-G [STR] | -B [STR] | -T [STR])  # Mutually exclusive args
        (-G [STR], --use_gisaid_data [STR]) <fasta_directory>
              Set flag if using single GISAID sequence file with single 
              GISAID metadata file. Indicate name of sequence file including 
              full path.
        (-B [STR], --use_genbank_data [STR]) <fasta_directory>
              Set flag if using single directory of GenBank fasta files. 
              Indicate full path of this directory.
        (-T [STR], --use_tax_tree [STR]) <fasta_directory>
              Set flag if using taxonomic tree structure with single genome 
              per species in lowest directory. Indicate path to root directory.
              **Not yet implemented**

    ('-r', '--resource_directory') <resource_directory>
        This is the directory which must contain certain necessary files
        and to which other resource files will be created.  Use the 
        absolute path.  See "evaluate_assay.py" in the evaluate_assay
        repo for details on which files need to be present.

    ('-R', '--results_directory') <results_directory>
        This is the directory to which results files will be created.
        Use the absolute path.

    ('-t', '--tnt_results_directory') <tnt_results_directory>
        This is the directory to which results files will be created when
        TNTBLAST is run.  Use the absolute path.
            
        
This script has the following optional arguments:
    
    (-A [STR], --assay_type [STR])
            Assay format to simulate
        
    (-e [INT], --min_primer_Tm [INT])
            Minimum primer Tm
        
    (-E [INT], --min_probe_Tm [INT])
            Minimum probe Tm
        
    (-l [INT], --max_amplicon_length [INT])
            Maximum amplicon length
        
    (-m [INT], --minimum_target_length [INT])
            Minimum length of target sequences. Doesn't filter GISAID
        
    [--best_match]
            If set, only save the best match, in Tm, between a query and target.

    [--tntblast_location [STR]]
            path to tntblast if system can't find it

@author: adanm
"""
import os
import sys
import time
import global_txt
import config_prep
import gbk_file_prep
import argparse as ap

''' Methods '''

def create_filename_namespace(args):
    filenames_ns = ap.Namespace()
    
    # Filenames
    config_filename = "variable_config.py"
    assay_filename = "assays.txt"
    gisaid_filename = "sequences.fasta"
    genbank_filename = "All_Seqs.fasta"
    if args.use_gisaid_data:
        sequence_filename = gisaid_filename
    else:
        sequence_filename = genbank_filename
    
    # Paths
    filenames_ns.resource_dir = args.resource_directory
    filenames_ns.results_dir = args.results_directory
    filenames_ns.fasta_dir = args.fasta_directory
    filenames_ns.tnt_dir = args.tnt_results_directory
        
    # Files
    filenames_ns.config_file = os.path.join(
        args.config_save_loc, config_filename)
    filenames_ns.assay_file = os.path.join(
        filenames_ns.resource_dir, assay_filename)
    filenames_ns.combined_fasta = os.path.join(
        filenames_ns.fasta_dir, sequence_filename)
    
    filenames_ns.db_stats_file = os.path.join(
        filenames_ns.results_dir, "db_stats.json")
    filenames_ns.db_totals_file = os.path.join(
        filenames_ns.results_dir, "db_totals.json")
    filenames_ns.accession_list_file = os.path.join(
        filenames_ns.results_dir, "accession_list.txt") 
    
    return filenames_ns


def prep_sequence_file(args, filenames_ns):
    
    if args.use_gisaid_data:
        print("...Using GISAID single download file.  No prep necessary")
        
    elif args.use_genbank_data:
        gbk_file_prep.prep_gbk(filenames_ns, args.minimum_target_length)
        
    elif args.use_tax_tree:
        print("Tax tree directory:", args.use_tax_tree)
        print("Taxonomic tree not yet implemented...")
    
    return


def get_arguments():
    parser = ap.ArgumentParser(prog=sys.argv[0].split('/')[-1])#, description=sys.argv[0].split('/')[-1])
    
    runtype_grp = parser.add_mutually_exclusive_group(required=True)
    runtype_grp.add_argument('-G', '--use_gisaid_data', action='store_true',
                        help="Set flag if using single GISAID sequence file" +
                        "with single GISAID metadata file.")
    runtype_grp.add_argument('-B', '--use_genbank_data', action='store_true',
                        help="Set flag if using single directory of GenBank " +
                        "fasta files.")
    runtype_grp.add_argument('-T', '--use_tax_tree', action='store_true',
                        help="Set flag if using taxonomic tree structure " +
                        "with single genome per species in lowest directory." +
                        "  Not yet implemented.")
    
    parser.add_argument('-A', '--assay_type', metavar='[STR]', type=str,
                        help="Assay format to simulate")
    parser.add_argument('-e', '--min_primer_Tm', metavar='[INT]', type=int,
                        help="Minimum primer Tm ")
    parser.add_argument('-E', '--min_probe_Tm', metavar='[INT]', type=int,
                        help="Minimum probe Tm ")
    parser.add_argument('-l', '--max_amplicon_length', metavar='[INT]', type=int,
                        help="Maximum amplicon length ")
    parser.add_argument('-m', '--minimum_target_length', metavar='[INT]', type=int,
                        default=0, help="Minimum length of target sequences. " +
                        "Doesn't filter GISAID")
    parser.add_argument('--best_match', action='store_true',
                        help="Only save the best match, in Tm, between a " +
                        "query and target")
    
    
    parser.add_argument('-c', '--config_save_loc', metavar='[STR]', type=str,
                        required=True, help="directory in which to save " +
                        "written config file")
    parser.add_argument('-r', '--resource_directory', metavar='[STR]', type=str,
                        required=True, help="resource files directory")
    parser.add_argument('-R', '--results_directory', metavar='[STR]', type=str,
                        required=True, help="results files directory")
    parser.add_argument('-f', '--fasta_directory', metavar='[STR]', type=str, 
                        required=True, help="location of fasta files")
    parser.add_argument('-t', '--tnt_results_directory', metavar='[STR]', type=str,
                        required=True, help="location of tnt results")
    parser.add_argument('--tntblast_location', metavar='[STR]', type=str,
                        help="path to tntblast if system can't find it")
    return parser.parse_args()



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    big_tic = time.perf_counter()

    ''' Get command line arguments '''
    args = get_arguments()
    # File Names Namespace
    filenames_ns = create_filename_namespace(args)
    
    
    
    ''' Write config file for tnt_wrapper to use '''
    print("\nWriting variable config file...")
    config_prep.write_config_file(args, filenames_ns, global_txt)
    
    
    ''' Prep target sequence file '''
    print("\nPrepping target sequence files...")
    prep_sequence_file(args, filenames_ns)
    
    
    
    




    big_toc = time.perf_counter()
    elapsed_time = big_toc - big_tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("\nFinished")
    print("Total Elapsed Time: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:12:38 2020

This script is a basic wrapper for ThermonucleotideBLAST (TNTBLAST).  
    For configuration it requires that the file variable_config.py is 
    present in the same directory.  For the time being, 
    variable_config.py is produced by the prep_for_assay_eval.py script 
    in the preprocess_assayeval repo (q.v. for details related to that 
    config file).  This will be modified to receive files from upstream 
    in the BioAI pipline at some point.


This script calls ThermonucleotideBLAST (TNTBLAST) which was created by 
    Jason Gans at Los Alamos National Laboratory with the BSD 3-Clause 
    License.  Please contact Jason at jgans@lanl.gov

    The reference for ThermonucleotideBLAST is: "Improved 
    assay-dependent searching of nucleic acid sequence databases." by 
    J. D. Gans and M. Wolinsky, Nucleic Acids Res. 2008 Jul;36(12):e74. 
    doi: 10.1093/nar/gkn301.

    Much thanks to Jason Gans for creating the algorithm that is the 
    core functionality for this package.

    To download and install TNTBLAST, please go to
    https://github.com/jgans/thermonucleotideBLAST


Dependencies and required files:
    TNTBLAST
    "variable_config.py"
    "assays.txt"

The following configuration file must be present in the same directory:
    
    "variable_config.py": Produced by and described in the 
            prep_for_assay_eval.py script in the 
            preprocess_assayeval repo
    
    
This file must be present in the resource directory (The resource
        directory is defined in the config):
    
    "assays.txt": A space-delimited (possibly tabs also?) text file that 
            contains a list of assays with corresponding oligo 
            sequences.  Each to a line.  Each line should contain the 
            following fields in order:

            "assay_name" "for_primer_seq" "rev_primer_seq" "probe_seq"
        Example line (newline escaped to fit comment size):
CDC-2019-nCoV_N1 GACCCCAAAATCAGCGAAAT TCTGGTTACTGCCAGTTGAATCTG \
    ACCCCGCATTACGTTTGGTGGACC

@author: adanm
"""
import os
import sys
import time
import argparse as ap
import variable_config as var_conf 


''' Methods '''

def create_filename_namespace(var_conf):
    filenames_ns = ap.Namespace()
    
    # Paths
    filenames_ns.resource_dir = var_conf.resource_directory
    filenames_ns.results_dir = var_conf.results_directory
    filenames_ns.tnt_dir = var_conf.tnt_results_directory
    filenames_ns.fasta_dir = var_conf.fasta_directory
    
    # Files
    filenames_ns.assay_file = var_conf.assay_file
    filenames_ns.combined_fasta = var_conf.sequence_file
    
    return filenames_ns


def prep_tnt_command(path_ns, var_ns):
    
    assay_file_arg = " -i {}".format(path_ns.assay_file)
    seq_file_arg = " -d {}".format(path_ns.combined_fasta)
    
    seq_file_prefix = path_ns.combined_fasta.split(".")[0]
    seq_acc = seq_file_prefix.split("/")[-1]
    tnt_results_file = "{}_results.out".format(seq_acc)
    tnt_results_path = os.path.join(path_ns.tnt_dir, tnt_results_file)
    output = " -o {}".format(tnt_results_path)
    
    runline_fields = []
    runline_fields.append(os.path.join(var_ns.tntblast_location, "tntblast") 
                            if var_ns.tntblast_location else "tntblast")
    runline_fields.append(" -A {}".format(var_ns.assay_type) 
                            if var_ns.assay_type else "")
    runline_fields.append(assay_file_arg)
    runline_fields.append(seq_file_arg)
    runline_fields.append(" -e {}".format(var_ns.min_primer_Tm) 
                            if var_ns.min_primer_Tm else "")
    runline_fields.append(" -E {}".format(var_ns.min_probe_Tm) 
                            if var_ns.min_probe_Tm else "")
    runline_fields.append(" -l {}".format(var_ns.max_amplicon_length) 
                            if var_ns.max_amplicon_length else "")
    runline_fields.append(output)
    runline_fields.append(" --best-match" if var_ns.best_match else "")

    tnt_runline = "".join(runline_fields)

    # return tnt_runline, tnt_results_path
    return tnt_runline


def run_tntblast_cmd(tnt_runline):
    os.system("{}".format(tnt_runline))


def run_tntblast(path_ns, var_ns):

    tnt_runline = prep_tnt_command(path_ns, var_ns)
    run_tntblast_cmd(tnt_runline)


    # tnt_results_list = []
    # tnt_runline, tnt_results = prep_tnt_command(path_ns, var_ns)
    # run_tntblast_cmd(tnt_runline)
    # tnt_results_list.append(tnt_results)
    
    # ''' ############### Output TNT Results List ################## '''
    # if args.use_gisaid_data:
    #     with open(tnt_results_list_file, 'w') as write_file:
    #         for line in tnt_results_list:
    #             print(line, file=write_file)
    return



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    big_tic = time.perf_counter()
    
    
    ''' Prep namespace '''
    path_ns = create_filename_namespace(var_conf)
    
    ''' Run TNTBLAST '''
    print("\nRunning tntblast against all sequences in combined file...")
    run_tntblast(path_ns, var_conf)


    big_toc = time.perf_counter()
    elapsed_time = big_toc - big_tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("\nFinished")
    print("Total Elapsed Time: {} hours, {} minutes, {} seconds".format(
            hours, minutes, seconds))


if __name__ == "__main__":
    main()

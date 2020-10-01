#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:58:09 2020

This script will download and prepare the necessary data for running the
    assay_monitor.py script.  Other scripts are called which access NCBI
    using NCBI's e-Utilities to download fasta files and metadata from 
    GenBank.  


Dependencies:
    ncbi_download.py


The script requires the following arguments:

        '-e', '--email' <email_address>'
            This is a requirement from ncbi for inclusion in a query

        '-r', '--resource_directory' <resource_directory>
            This is the directory to which certain resource files will be
            created.
            Use the format:  /XXX/XXX/XXX/

        '-R', '--results_directory' <results_directory>
            This is the directory to which certain results files will be
            downloaded or created.
            Use the format:  /XXX/XXX/XXX/

        '-f', '--fasta_directory' <download_directory>
            This is the directory to which fna sequence files will be
            downloaded.
            Use the format:  /XXX/XXX/XXX/

        '-m', '--mindate' <early_date>
            Put the earliest publication date or lower range of date in format:
            "YYYY/MM/DD"

        '-M', '--maxdate' <late_date>
            Put the latest publication date or upper range of date in format:
            "YYYY/MM/DD"

        '-s', '--seqlengthrange' <seq_length_range>
            Input the range of sequence lengths in the format:  29000:30500.
            For no limit in either direction, put: 0:99999999

@author: adanm
"""

import os
import sys
import time
import argparse as ap
import ncbi_download





def get_arguments():
    parser = ap.ArgumentParser(prog=sys.argv[0].split('/')[-1], description=sys.argv[0].split('/')[-1])
    parser.add_argument('-e', '--email', metavar='[STR]', nargs=1, type=str,
                        required=True, help="return email from ncbi")
    parser.add_argument('-m', '--mindate', metavar='[STR]', nargs=1, type=str,
                        required=True, help="minimum publication date")
    parser.add_argument('-M', '--maxdate', metavar='[STR]', nargs=1, type=str,
                        required=True, help="maximum publication date")
    parser.add_argument('-s', '--seqlengthrange', metavar='[STR]', nargs=1, type=str,
                        required=True, help="range of sequence lengths")
    parser.add_argument('-r', '--resource_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="resource files directory")
    parser.add_argument('-R', '--results_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="results files directory")
    parser.add_argument('-f', '--fasta_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="location of fasta files")

    return parser.parse_args()



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    big_tic = time.perf_counter()

    ''' Get command line arguments '''
    args = get_arguments()


    tool = sys.argv[0].split('/')[-1]
    email = args.email[0]

    min_date = args.mindate[0]
    max_date = args.maxdate[0]
    seq_lengths = args.seqlengthrange[0]


    resource_dir = args.resource_directory[0]
    results_dir = args.results_directory[0]
    fasta_dir = args.fasta_directory[0]


    # Query
    # target_query = 'txid2697049+AND+{}[SLEN]'.format(seq_lengths)  # alternative query
    target_query = 'SARS-cov-2[orgn]+AND+{}[SLEN]'.format(seq_lengths)


    # File Names
    esearch_query_file = os.path.join(resource_dir, 'esearch_results.xml')
    efetch_results_file = os.path.join(results_dir, 'efetch_results.txt')


    ''' ESearch '''
    print("5 second pause...")
    time.sleep(5)
    ncbi_download.run_esearch(target_query, esearch_query_file, email, tool, min_date, max_date)


    ''' Get id list from esearch file '''
    id_list = ncbi_download.get_id_list(esearch_query_file)
    print("Length of search list:", len(id_list), "\n")

    ''' Get indices to split list '''
    indices = ncbi_download.split_indices(id_list)


    ''' EFetch '''
    for index in indices:
        print("\nGetting sequences {} to {}...".format(index[0], index[1]))
        id_list_string = ncbi_download.split_id_list(id_list, index)


        print("5 second pause...")
        time.sleep(5)
        ncbi_download.run_efetch(id_list_string, fasta_dir, efetch_results_file, email, tool)

        file_list = ncbi_download.read_efetch_write_fna(efetch_results_file, fasta_dir)
        for each in file_list:
                print(each)


    big_toc = time.perf_counter()
    elapsed_time = big_toc - big_tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60 - hours * 3600
    print("\nFinished")
    print("Total Elapsed Time: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))


if __name__ == "__main__":
    main()

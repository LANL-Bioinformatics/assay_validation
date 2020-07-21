#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:58:09 2020

This script will download and prepare the necessary data for running the
    assay_monitor.py script.  Other scripts are called which log into the
    GISAID database to download SARS-CoV-2 genomes as individual fasta (.fna)
    files.  NCBI is also accessed using NCBI's e-Utilities to download fastas
    from GenBank.  Metadata for both GISAID and GenBank are downloaded from
    these sources.


This script calls gisaid_EpiCoV_downloader.py which was created by Paul Li at
    LANL in 2020 with the GPL 1.0.0 license.  Contact Paul at po-e@lanl.gov

This script requires that you have access to GISAID to be able to log in and
    download sequences.

We gratefully acknowledge the team at GISAID for creating the the COVID-19 outbreak global database and resources, and the many authors from the originating and submitting laboratories of for the SARS-CoV-2 sequence data, which provides the foundation for the analyses tools provided here. The original data are available from https://www.gisaid.org (registration required). A full acknowledgements of all originating and submitting laboratories can be downloaded.


Dependencies:
    gisaid_split.py
    ncbi_download.py
    gisaid_EpiCoV_downloader.py


The script requires the following arguments:

        '-u', '--gisaid_username' <username>
            This is your username to log into gisaid.org

        '-p', '--gisaid_passwd' <password>
            This is your password for your GISAID login

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

        '-t', '--tnt_results_directory' <tnt_results_directory>
            This is the directory to which results files will be created when
            TNTBLAST is run.
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

import re
import os
import sys
import time
import test
import gisaid_split
import subprocess

import argparse as ap
import ncbi_download
import gisaid_EpiCoV_downloader as gis_downldr





def get_arguments():
    parser = ap.ArgumentParser(prog=sys.argv[0].split('/')[-1], description=sys.argv[0].split('/')[-1])
    parser.add_argument('-u', '--gisaid_username', metavar='[STR]', nargs=1, type=str,
                        required=True, help="login name for gisaid")
    parser.add_argument('-p', '--gisaid_passwd', metavar='[STR]', nargs=1, type=str,
                        required=True, help="password for gisaid login")
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
    parser.add_argument('-t', '--tnt_results_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="location of tnt results")
    return parser.parse_args()



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    big_tic = time.perf_counter()

    ''' Get command line arguments '''
    args = get_arguments()

    gname = args.gisaid_username[0]
    gpass = args.gisaid_passwd[0]

    tool = sys.argv[0].split('/')[-1]
    email = args.email[0]

    min_date = args.mindate[0]
    max_date = args.maxdate[0]
    seq_lengths = args.seqlengthrange[0]

    resource_dir = args.resource_directory[0]
    results_dir = args.results_directory[0]
    fasta_dir = args.fasta_directory[0]
    tnt_res_dir = args.tnt_results_directory[0]


    # Query
    # target_query = 'txid2697049+AND+{}[SLEN]'.format(seq_lengths)  # alternative query
    target_query = 'SARS-cov-2[orgn]+AND+{}[SLEN]'.format(seq_lengths)


    # File Names
    sequences_file = os.path.join(fasta_dir, "sequences.fasta")
    sequences_download_zip = os.path.join(fasta_dir, "sequences*.fasta.gz")
    metadata_file = os.path.join(resource_dir, "metadata.tsv")
    metadata_download_zip = os.path.join(fasta_dir, "metadata*.tsv.gz")

    metafile = os.path.join(resource_dir, "metadata.tsv")
    gisaid_fasta_file = os.path.join(fasta_dir, "sequences.fasta")

    esearch_query_file = os.path.join(resource_dir, 'esearch_results.xml')
    efetch_results_file = os.path.join(results_dir, 'efetch_results.txt')



    ''' ################ Download GISAID Seqs and Metadata ################# '''

    gis_downldr.download_gisaid_EpiCoV(uname=gname, upass=gpass,
        normal=False, wd=fasta_dir,
        loc=None,
        cs=None,
        ce=None,
        ss=None,
        se=None,
        cg=False,
        hc=False,
        le=False,
        to=90,
        rt=5,
        iv=3,
        meta_dl=False )


    if os.path.exists(sequences_file):
        os.remove(sequences_file)

    if os.path.exists(metadata_file):
        os.remove(metadata_file)

    runline = "gunzip {}".format(sequences_download_zip)
    os.system("{}".format(runline))

    sequences_download = [ os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if re.match(r"sequences.*.fasta", f) ][-1]

    os.rename(sequences_download, sequences_file)

    runline = "gunzip {}".format(metadata_download_zip)
    os.system("{}".format(runline))


    metadata_download = [ os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if re.match(r"metadata.*.tsv", f) ][-1]

    os.rename(metadata_download, metadata_file)



    ''' ########################## GISAID Split ############################ '''


    ''' Split gisaid fasta into multiple fna files in target directory '''
    if os.path.exists(gisaid_fasta_file):
        print("\nSplitting GISAID sequence file...")
        gisaid_split.split_gisaid_file_with_meta(gisaid_fasta_file, fasta_dir, metafile)
    else:
        print("\nNo Gisaid Fasta File !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")



    ''' ######################### GenBank Search & Fetch ##################### '''

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
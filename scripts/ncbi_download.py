#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 13:42:33 2020

@author: adanm
"""
import os
import sys
import time
import requests
import argparse as ap
import xml.etree.ElementTree as ET


''' 'Global' Variables '''
# Data download switches
post_esearch = True
post_efetch = True
write_files = True

# Search Elements
base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
db = 'nuccore'
use_history = 'n'
return_type = 'null' # text ASN.1
return_mode = 'text' # text ASN.1
return_type = 'native' # Full record in XML
return_type = 'fasta' # TinySeq XML
return_type = 'gb' # GBSeq XML
return_type = 'gbc' # INSDSeq XML
return_mode = 'xml'
max_returned_recs = 3000
split_size = 200
date_type = 'pdat'



''' Methods '''

def run_esearch(query, search_results_file, email, tool, min_date, max_date):

    # assemble the esearch URL
    url = base_url + "esearch.fcgi?db={}&term={}&usehistory={}&retmax={}".format(
            db, query, use_history, max_returned_recs)
    url = url + "&mindate={}&maxdate={}&datetype={}".format(min_date, max_date, date_type)
    url = url + "&tool={}&email={}".format(tool, email)

    #post the esearch URL
    print("\nPosting following search:")
    print(url, "\n")
    search_results = requests.get(url)

    # Write to file
    with open(search_results_file, 'w') as file_handle:
        file_handle.write(search_results.text)

    print("eSearch completed")
    return


def get_id_list(search_results_file):
    # parse results file for IDs
    tree = ET.parse(search_results_file)
    root = tree.getroot()

    # Make ID list
    id_list = []
    for id_num in root.iter('Id'):
        id_list.append(id_num.text)

    return id_list


def split_indices(id_list):
    ''' Get indices for splitting id list '''
    list_length = len(id_list)

    split_idx_list = []
    if list_length >= split_size:
        splits = int( list_length / split_size )
        for i in range(splits):
            split_idx_list.append([ (i)*split_size, (i + 1)*split_size ])
        if (i + 1)*split_size != list_length:
            split_idx_list.append([(i + 1)*split_size, list_length])
    else:
        split_idx_list.append([0, list_length])

    return split_idx_list


def split_id_list(id_list, index_pair):

    # Split id_list using index pair and return list as a string
    id_list = id_list[index_pair[0]:index_pair[1]]
    print("Length of split search list:", len(id_list), "\n")
    id_list_string = ','.join(id_list)

    return id_list_string


def run_efetch(id_list_string, fna_download_path, data_file, email, tool):

    #assemble the efetch URL
    url = base_url + "efetch.fcgi?db={}&id={}".format(db, id_list_string)
    url += "&rettype={}&retmode={}".format(return_type, return_mode)
    url = url + "&tool={}&email={}".format(tool, email)


    #post the efetch URL
    print("Posting following fetch from NCBI:")
    print(url, "\n")
    data = requests.get(url)

    # Write fetched data to file
    with open(data_file, 'w') as file_handle:
        file_handle.write(data.text)


def read_efetch_write_fna(data_file, sequence_dir):

    # Parse INSD xml file
    tree = ET.parse(data_file)
    root = tree.getroot()

    ''' Get Metadata and sequence for a single sequence from INSD xml file '''
    file_list = []
    # Initial Metadata
    for node in root.iter("INSDSeq"):
        locus = node.find("INSDSeq_locus").text

        seq_length =  node.find("INSDSeq_length").text
        create_date = node.find("INSDSeq_create-date").text
        accession_number = node.find("INSDSeq_primary-accession").text
        try:
            comment = node.find("INSDSeq_comment").text
        except:
            print(locus, "has no comment #####################################################")
            comment = ""

        # Feature Table metadata
        feature_table = node.find("INSDSeq_feature-table")
        # Interval
        for feature_node in feature_table.iter("INSDFeature"):
            feat_key = feature_node.find("INSDFeature_key").text
            if feat_key == "source":
                interval_node = feature_node.find("INSDFeature_intervals")
                next_node = interval_node.find("INSDInterval")
                int_start = next_node.find("INSDInterval_from").text
                int_end = next_node.find("INSDInterval_to").text
        # Other Metadata
        for qualifier_node in feature_table.iter("INSDQualifier"):
            qual_name = qualifier_node.find("INSDQualifier_name").text
            try:
                qual_value = qualifier_node.find("INSDQualifier_value").text
            except:
                # print(locus, qual_name, "has no information, skipping...")
                pass
            if qual_name == "isolate":
                isolate = qual_value
            if qual_name == "country":
                country = qual_value
            if qual_name == "collection_date":
                col_date = qual_value

        # Sequence
        try:
            sequence_text = node.find("INSDSeq_sequence").text
        except:
            sequence_text = ""
            print("No sequence in ", locus, "#######################################################")


        # Create header using metadata
        header_text = ">{}|{}|{}|{}|{}|{}|{}|{}|{}".format(
            isolate, locus, col_date, create_date, country, seq_length,
            int_start, int_end, comment)

        # Write to fna file
        output_fna_file = os.path.join(sequence_dir, "{}.fna".format(accession_number))
        write_fna(output_fna_file, header_text, sequence_text.upper())
        file_list.append(output_fna_file)

    return file_list


def write_fna(file_name, header, sequence):
    with open(file_name, 'w') as file_handle:
        file_handle.write(header + "\n")
        file_handle.write(sequence +"\n")


def get_arguments():
    parser = ap.ArgumentParser(prog='ncbi_download.py',
                               description="""NCBI download script""")
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
                        required=True, help="fasta download directory")
    return parser.parse_args()


def main():
    ''' Get command line argument '''
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


    ''' ######################### GenBank Search & Fetch ##################### '''

    ''' ESearch '''
    if post_esearch:
        print("5 second pause...")
        time.sleep(5)
        run_esearch(target_query, esearch_query_file, email, tool, min_date, max_date)


    ''' Get id list from esearch file '''
    id_list = get_id_list(esearch_query_file)
    print("Length of search list:", len(id_list), "\n")

    ''' Get indices to split list '''
    indices = split_indices(id_list)


    ''' EFetch '''
    for index in indices:
        print("\nGetting sequences {} to {}...".format(index[0], index[1]))
        id_list_string = split_id_list(id_list, index)

        if post_efetch:
            print("5 second pause...")
            time.sleep(5)
            run_efetch(id_list_string, fasta_dir, efetch_results_file, email, tool)

        if write_files:
            file_list = read_efetch_write_fna(efetch_results_file, fasta_dir)
            for each in file_list:
                    print(each)

    print()
    print("Finished")


if __name__ == "__main__":
    main()

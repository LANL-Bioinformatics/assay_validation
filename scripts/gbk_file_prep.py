#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 17:22:40 2021

@author: adanm
"""
import re
import os
import db_stats
import seq_db_filter_funx as filt_db


''' Methods '''

def generate_fna_lists(fasta_dir):

    ''' Get list of all genbank files '''
    indiv_file_list = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir)
                 if re.match(r'.*\.fna', f) ]    
    indiv_file_list.sort()

    return indiv_file_list


def prep_gbk(filenames_ns, minimum_file_length_considered):
    
    fasta_dir = filenames_ns.fasta_dir
    combined_fasta = filenames_ns.combined_fasta
    db_stats_file = filenames_ns.db_stats_file
    db_totals_file = filenames_ns.db_totals_file
    accession_list_file = filenames_ns.accession_list_file
    
    
    ''' ####################### Generate Databases ######################### '''

    genbank_file_list = generate_fna_lists(fasta_dir)

    ''' ######################## Filter Databases ########################## '''

    print("Filtering GenBank...")
    file_list, genbank_tallies, total_list_len = filt_db.filter_len_et_orgs(
        genbank_file_list, minimum_file_length_considered)

    ''' ##################### Combine Sequence Files ####################### '''

    ''' Combine fna sequence files into single fasta file '''
    print("Combining target fna files from database into one fasta...")
    accession_list = []
    with open(combined_fasta, 'w') as write_combine:
        for fna_file in file_list:
            accession = fna_file.split('/')[-1].split('.')[0]

            with open(fna_file, 'r') as read_fna:
                for line in read_fna.readlines():
                    print(line.strip(), file=write_combine)
            accession_list.append(accession)
            
    ''' ##################### Output Accession List ######################## '''
    with open(accession_list_file, 'w') as write_file:
        for line in accession_list:
            print(line, file=write_file)


    ''' ##################### Output DB Stats Files ######################## '''
    db_stats.generate_db_stats(genbank_tallies, total_list_len, db_stats_file, db_totals_file)
        
    db_stats.generate_db_stats(genbank_tallies, total_list_len, db_stats_file, db_totals_file)
    
    return



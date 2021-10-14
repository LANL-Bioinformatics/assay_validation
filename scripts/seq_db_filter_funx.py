#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 16:18:46 2020

@author: adanm
"""

def check_fna_file(fna_file, minimum_file_length_considered):
    with open(fna_file, 'r') as file_handle:
        line1 = file_handle.readline().strip()
        line2 = file_handle.readline().strip()
        remaining_length = len(file_handle.readlines())
        if remaining_length > 0:
            print("length of remaining lines:", remaining_length, fna_file)
    if "bat" in line1:
        return False
    if "pangolin" in line1:
        return False
    if len(line2) < minimum_file_length_considered:
        return False
    return True


def tally_removed(fna_file, minimum_file_length_considered, tally_list):
    with open(fna_file, 'r') as file_handle:
        line1 = file_handle.readline().strip()
        line2 = file_handle.readline().strip()
    if "bat" in line1:
        tally_list[0] += 1
    elif "pangolin" in line1:
        tally_list[1] += 1
    elif len(line2) < minimum_file_length_considered:
        tally_list[2] += 1
    return tally_list



def filter_len_et_orgs(file_list, min_seq_len):
    bat_pang_short_tally = [ 0, 0, 0]
    removed = 0
    kept = 0
    list_to_return = file_list.copy()
    list_before = len(list_to_return)
    for fna_file in file_list:
        accession = fna_file.split('/')[-1].split('.')[0]
        if not check_fna_file(fna_file, min_seq_len):
            print("Omitted file {} (Failed check) ".format(accession) + \
                  "#####################################################")
            bat_pang_short_tally = tally_removed(fna_file, min_seq_len, bat_pang_short_tally)
            list_to_return.remove(fna_file)
            removed += 1
        else:
            # print("{} checked".format(accession))
            kept += 1
    list_after = len(list_to_return)

    tallies = [ list_before, removed, bat_pang_short_tally[0], 
               bat_pang_short_tally[1], bat_pang_short_tally[2], 
               kept, list_after ]
    
    total_list_len = len(list_to_return)

    return list_to_return, tallies, total_list_len


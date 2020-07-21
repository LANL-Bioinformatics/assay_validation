#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 16:18:46 2020

@author: adanm
"""
Remove_GISAID_Duplicates = False # else remove GenBank duplicates

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


def filter_overlap(genbank_file_list, gisaid_file_list):
    genbank_start_number = len(genbank_file_list)
    gisaid_start_number = len(gisaid_file_list)

    ''' Get isolate metadata from genbank files '''
    isolate_field_list = []
    for gen_file in genbank_file_list:
        with open(gen_file, 'r') as read_genfile:
            header = read_genfile.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            isolate_fields = isolate.split('/')
            isolate_field_list.append(isolate_fields)

    ''' From gisaid file list: Remove overlap with genbank '''
    filtered_gisaid_list = gisaid_file_list.copy()

    overlap = 0
    exception = 0
    for gis_file in gisaid_file_list:
        with open(gis_file, 'r') as read_gisfile:
            header = read_gisfile.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            isolate_fields = isolate.split('/')
            if len(isolate_fields) >=4:
                for gis_isolate in isolate_field_list:
                    if isolate_fields[1] in gis_isolate and isolate_fields[2] in gis_isolate:
                        try:
                            print("Removing", gis_file)
                            filtered_gisaid_list.remove(gis_file)
                            overlap += 1
                        except:
                            exception += 1
                            pass
            else:
                for gis_isolate in isolate_field_list:
                   if isolate_fields[1] in gis_isolate:
                        try:
                            print("Removing", gis_file)
                            filtered_gisaid_list.remove(gis_file)
                            overlap += 1
                        except:
                            exception += 1
                            pass
    gisaid_remove_overlap = len(filtered_gisaid_list)

    tallies = [ genbank_start_number, gisaid_start_number, exception, overlap, gisaid_remove_overlap]


    return genbank_file_list, filtered_gisaid_list, tallies


def filter_overlap_remove_genbank(genbank_file_list, gisaid_file_list):
    genbank_start_number = len(genbank_file_list)
    gisaid_start_number = len(gisaid_file_list)

    ''' Get isolate metadata from genbank files '''
    isolate_field_list = []
    for gis_file in gisaid_file_list:
        with open(gis_file, 'r') as read_genfile:
            header = read_genfile.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            isolate_fields = isolate.split('/')
            isolate_field_list.append(isolate_fields)


    ''' From gisaid file list: Remove overlap with genbank '''
    filtered_genbank_list = genbank_file_list.copy()

    overlap = 0
    exception = 0
    for gen_file in genbank_file_list:
        with open(gen_file, 'r') as read_genfile:
            header = read_genfile.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            isolate_fields = isolate.split('/')
            if len(isolate_fields) >=4:
                for gis_isolate in isolate_field_list:
                    if isolate_fields[2] in gis_isolate and isolate_fields[3] in gis_isolate:
                        try:
                            print("Removing", gen_file)
                            filtered_genbank_list.remove(gen_file)
                            overlap += 1
                        except:
                            exception += 1
                            pass
            elif len(isolate_fields) ==3:
                for gis_isolate in isolate_field_list:
                    # print("\n", isolate_fields)
                    # print(gis_isolate, isolate_fields[1])
                    if isolate_fields[1] in gis_isolate:
                        try:
                            print("Removing", gen_file)
                            filtered_genbank_list.remove(gen_file)
                            overlap += 1
                        except:
                            exception += 1
                            pass
            elif len(isolate_fields) ==1:
                for gis_isolate in isolate_field_list:
                    # print("\n", isolate_fields)
                    # print(gis_isolate, isolate_fields[0])
                    if isolate_fields[0] in gis_isolate:
                        try:
                            print("Removing", gen_file)
                            filtered_genbank_list.remove(gen_file)
                            overlap += 1
                        except:
                            exception += 1
                            pass
    genbank_remove_overlap = len(filtered_genbank_list)

    tallies = [ genbank_start_number, gisaid_start_number, exception, overlap, genbank_remove_overlap]


    return filtered_genbank_list, gisaid_file_list, tallies


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

    tallies = [ list_before, removed, bat_pang_short_tally[0], bat_pang_short_tally[1], bat_pang_short_tally[2], kept, list_after ]

    return list_to_return, tallies


def filter_dbs(genbank_file_list, gisaid_file_list, min_seq_len):

    print("\nFiltering Overlap...")
    genbank_file_list, gisaid_file_list, overlap_tallies = filter_overlap(genbank_file_list, gisaid_file_list)

    print("\nFiltering GenBank...")
    genbank_file_list, genbank_tallies = filter_len_et_orgs(genbank_file_list, min_seq_len)

    print("\nFiltering GISAID...")
    gisaid_file_list, gisaid_tallies = filter_len_et_orgs(gisaid_file_list, min_seq_len)

    ''' Combine GenBank and GISAID '''
    file_list = gisaid_file_list + genbank_file_list
    total_list_len = len(file_list)

    return file_list, overlap_tallies, genbank_tallies, gisaid_tallies, total_list_len



def filter_dbs_overlap_last(genbank_file_list, gisaid_file_list, min_seq_len):


    print("\nFiltering GenBank...")
    genbank_file_list, genbank_tallies = filter_len_et_orgs(genbank_file_list, min_seq_len)

    print("\nFiltering GISAID...")
    gisaid_file_list, gisaid_tallies = filter_len_et_orgs(gisaid_file_list, min_seq_len)

    print("\nFiltering Overlap...")
    if Remove_GISAID_Duplicates:
        genbank_file_list, gisaid_file_list, overlap_tallies = filter_overlap(
            genbank_file_list, gisaid_file_list)
    else:
        genbank_file_list, gisaid_file_list, overlap_tallies = filter_overlap_remove_genbank(
            genbank_file_list, gisaid_file_list)

    ''' Combine GenBank and GISAID '''
    file_list = gisaid_file_list + genbank_file_list
    total_list_len = len(file_list)

    return file_list, overlap_tallies, genbank_tallies, gisaid_tallies, total_list_len



def get_overlap_list(genbank_file_list, gisaid_file_list):
    genbank_start_number = len(genbank_file_list)
    gisaid_start_number = len(gisaid_file_list)

    ''' Get isolate metadata from genbank files '''
    isolate_field_dict = {}
    isolate_field_meta_dict = {}
    for gen_file in genbank_file_list:
        with open(gen_file, 'r') as read_genfile:
            header = read_genfile.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            isolate_fields = isolate.split('/')
            isolate_field_dict[header_fields[1]] = isolate_fields
            isolate_field_meta_dict[header_fields[1]] = header_fields
            # isolate_field_list.append(isolate_dict)

    ''' From gisaid file list: Remove overlap with genbank '''
    filtered_gisaid_list = gisaid_file_list.copy()

    overlap = 0
    exception = 0
    overlap_list = []
    for gis_file in gisaid_file_list:
        with open(gis_file, 'r') as read_gisfile:
            header = read_gisfile.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            accession = header_fields[1]
            isolate_fields = isolate.split('/')
            if len(isolate_fields) >=4:
                for gis_accession, gis_isolate in isolate_field_dict.items():
                    if isolate_fields[1] in gis_isolate and isolate_fields[2] in gis_isolate:
                        try:
                            table_line = [accession, isolate, gis_accession, "/".join(gis_isolate)]
                            if len(header_fields) < 18:
                                while len(header_fields) < 18:
                                    header_fields.append("")
                            table_line.extend(header_fields[2:])
                            table_line.extend(isolate_field_meta_dict[gis_accession][2:])
                            overlap_list.append(table_line)
                            overlap += 1
                        except:
                            exception += 1
                            pass
            else:
                for gis_accession, gis_isolate in isolate_field_dict.items():
                   if isolate_fields[1] in gis_isolate:
                        try:
                            table_line = [accession, isolate, gis_accession, "/".join(gis_isolate)]
                            if len(header_fields) < 18:
                                while len(header_fields) < 18:
                                    header_fields.append("")
                            table_line.extend(header_fields[2:])
                            table_line.extend(isolate_field_meta_dict[gis_accession][2:])
                            overlap_list.append(table_line)
                            overlap += 1
                        except:
                            exception += 1
                            pass
    gisaid_remove_overlap = len(filtered_gisaid_list)

    tallies = [ genbank_start_number, gisaid_start_number, exception, overlap, gisaid_remove_overlap]


    return overlap_list

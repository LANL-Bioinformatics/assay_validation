#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 13:42:33 2020

@author: adanm
"""
import os
import requests
import xml.etree.ElementTree as ET


''' 'Global' Variables '''
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


def read_efetch_write_fna_legacy(data_file, sequence_dir):

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
        isolate = ""
        country = ""
        col_date = ""
        comment = ""

        comment_node = node.find("INSDSeq_comment")
        if comment_node is not None:
            comment = comment_node.text

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
            qual_value_node = qualifier_node.find("INSDQualifier_value")
            if qual_value_node is not None:
                qual_value = qual_value_node.text
            else:
                print(locus, qual_name, "has no information, skipping...")

            if qual_name == "isolate":
                isolate = qual_value
            if qual_name == "country":
                country = qual_value
            if qual_name == "collection_date":
                col_date = qual_value

        # Sequence
        sequence_text_node = node.find("INSDSeq_sequence")
        if sequence_text_node is not None:
            sequence_text = sequence_text_node.text
        else:
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


def read_efetch_write_fna(data_file, sequence_dir):

    # Parse INSD xml file
    tree = ET.parse(data_file)
    root = tree.getroot()

    ''' Get Metadata and sequence for a single sequence from INSD xml file '''
    file_list = []
    for node in root.iter("INSDSeq"):

        # Initialize
        segment = ""
        isolate = ""
        host = ""
        txid = ""
        country = ""
        division = ""
        country_field = ""
        col_date = ""
        comment = ""
        title = ""
        submitting_lab = ""

        # Get Basic Metadata
        locus = node.find("INSDSeq_locus").text
        seq_length =  node.find("INSDSeq_length").text
        create_date = node.find("INSDSeq_create-date").text
        common_name = node.find("INSDSeq_definition").text
        accession_number = node.find("INSDSeq_primary-accession").text
        version = node.find("INSDSeq_accession-version").text
        taxonomy = node.find("INSDSeq_taxonomy").text
        comment_node = node.find("INSDSeq_comment")
        if comment_node is not None:
            comment = comment_node.text


        # Reference metadata
        references = node.find("INSDSeq_references")
        for reference in references.iter("INSDReference"):
            ref_ref = reference.find("INSDReference_reference")
            title_node = reference.find("INSDReference_title")
            j_node = reference.find("INSDReference_journal")
            title_text = title_node.text
            j_title = j_node.text

            if title_text == "Direct Submission":
                submitting_lab = j_title
            else:
                title = title_text

            author_list = []
            for author in reference.iter("INSDAuthor"):
                author_list.append(author.text)
            authors = ",".join(author_list)


        # Feature Table metadata
        feature_table = node.find("INSDSeq_feature-table")
        for feature_node in feature_table.iter("INSDFeature"):
            feat_key = feature_node.find("INSDFeature_key").text
            if feat_key == "source":
                # Interval Metadata
                interval_node = feature_node.find("INSDFeature_intervals")
                next_node = interval_node.find("INSDInterval")
                int_start = next_node.find("INSDInterval_from").text
                int_end = next_node.find("INSDInterval_to").text
                # Other Metadata
                for qualifier_node in feature_node.iter("INSDQualifier"):
                    qual_name = qualifier_node.find("INSDQualifier_name").text
                    qual_value_node = qualifier_node.find("INSDQualifier_value")
                    if qual_value_node is not None:
                        qual_value = qual_value_node.text
                    else:
                        qual_value = ""

                    if qual_name == "mol_type":
                        segment = qual_value
                    if qual_name == "isolate":
                        isolate = qual_value
                    if qual_name == "host":
                        host = qual_value
                    if qual_name == "db_xref":
                        txid = qual_value.split(':')[-1]
                    if qual_name == "country":
                        country_field = qual_value
                        if ": " in country_field:
                            country_fields = country_field.split(": ")
                            country = country_fields[0]
                            division = country_fields[1]
                        elif ":" in country_field:
                            country_fields = country_field.split(":")
                            country = country_fields[0]
                            division = country_fields[1]
                        else:
                            country = country_field
                    if qual_name == "collection_date":
                        col_date = qual_value

        # Sequence
        sequence_text_node = node.find("INSDSeq_sequence")
        if sequence_text_node is not None:
            sequence_text = sequence_text_node.text
        else:
            sequence_text = ""
            print("No sequence in ", locus, "#######################################################")


        # Initialize empty fields.  Hopefully find entries for these fields some day
        region = ""
        location = ""
        region_exposure = ""
        country_exposure = ""
        division_exposure = ""
        age = ""
        sex = ""
        originating_lab = ""


        # Create header using metadata
        header_text = ">{}|{}|{}|{}|{}|{}|{}|{}|{}/{}/{}|{}/{}/{}/{}|{}|{}|{}|{}|{}|{}|{}|{}|{}".format(
            isolate, locus, col_date, create_date, country_field, seq_length, int_start,
            int_end, region, country, division, location, region_exposure,
            country_exposure, division_exposure, segment, host, age, sex,
            originating_lab, submitting_lab, authors, title, comment)

        # Write to fna file
        output_fna_file = os.path.join(sequence_dir, "{}.fna".format(accession_number))
        write_fna(output_fna_file, header_text, sequence_text.upper())
        file_list.append(output_fna_file)

    return file_list


def write_fna(file_name, header, sequence):
    with open(file_name, 'w') as file_handle:
        file_handle.write(header + "\n")
        file_handle.write(sequence +"\n")


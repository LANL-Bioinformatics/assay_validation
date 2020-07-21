#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:36:54 2020

@author: adanm
"""

import re
import os
import json
import argparse as ap
from datetime import datetime

''' 'Global' Variables '''
# Dictionary Labels
assay = "Assay"
assay_name = "Name"
FP = "Forward Primer"
RP = "Reverse Primer"
Pr = "Probe"
gp = "gaps"
mm = "mismatches"


''' Methods '''

def generate_table(in_file, meta_filename, sum_data_file, table_dict, short_assay_list):

    ''' Get data list from table dictionary '''
    data_list = table_dict["data"]

    ''' Get assay monitor data from json file '''
    with open(in_file) as json_file:
        full_dict = json.load(json_file)

    ''' loop through all assays in assay monitor dictionary '''
    for idx, assay_dict in enumerate(full_dict[assay]):
        this_assay = assay_dict[assay_name]

        if this_assay in short_assay_list:

            ''' initialize tallies '''
            succeed_tally = 0
            succeed_one_two = 0
            one_mm = 0
            two_mm = 0
            two_plus = 0
            three_mm = 0
            three_plus_fail = 0
            four_mm = 0
            five_mm = 0
            six_mm = 0
            seven_mm = 0
            eight_plus_mm = 0
            eight_plus_fail = 0
            fail_tally = 0
            pos_elims = 0
            neg_elims = 0

            ''' Get mismatch tallies from true positive list '''
            for seq_entry in assay_dict["True Positives"]:
                try:
                    mismatches = 0
                    mismatches_fp = int(seq_entry["Values"]["Forward Primer"]['mismatches'])
                    if 'mismatches' in seq_entry["Values"]["Probe"].keys():
                        mismatches_pr = int(seq_entry["Values"]["Probe"]['mismatches'])
                    else:
                        mismatches_pr = 0
                    mismatches_rp = int(seq_entry["Values"]["Reverse Primer"]['mismatches'])
                    mismatches = max(mismatches_fp, mismatches_pr, mismatches_rp)
                    if mismatches == 0:
                        succeed_tally += 1
                    elif mismatches == 1:
                        one_mm += 1
                    elif mismatches == 2:
                        two_mm += 1
                    elif mismatches == 3:
                        three_mm += 1
                    elif mismatches == 4:
                        four_mm += 1
                    elif mismatches == 5:
                        five_mm += 1
                    elif mismatches == 6:
                        six_mm += 1
                    elif mismatches == 7:
                        seven_mm += 1
                    elif mismatches >= 8:
                        eight_plus_mm += 1
                        eight_plus_fail += 1
                    if mismatches >= 2:
                        two_plus += 1
                    if mismatches >= 3:
                        three_plus_fail += 1
                    if mismatches < 3:
                        succeed_one_two += 1
                except:
                    pos_elims += 1

            ''' Get tallies of false negatives '''
            for seq_entry in assay_dict["False Negatives"]:
                try:
                    three_plus_fail += 1
                    eight_plus_fail += 1
                    fail_tally += 1
                except:
                    neg_elims += 1

            ''' Calculate Recall '''
            try:
                TPR = str(succeed_one_two/(succeed_one_two + three_plus_fail))
            except:
                print("this_assay", this_assay)
                print("succeed_one_two", succeed_one_two)
                print("three_plus_fail", three_plus_fail)


            ''' Prepare dictionary entry '''
            dict_entry = {}
            dict_entry["name"] = this_assay
            dict_entry["recall"] = TPR
            dict_entry["perfect_match"] = succeed_tally
            dict_entry["match_1_2"] = succeed_one_two
            dict_entry["1_mm"] = one_mm
            dict_entry["2_mm"] = two_mm
            dict_entry["2_mm_p"] = two_plus
            dict_entry["3_mm"] = three_mm
            dict_entry["3_mm_p_fail"] = three_plus_fail
            dict_entry["4_mm"] = four_mm
            dict_entry["5_mm"] = five_mm
            dict_entry["6_mm"] = six_mm
            dict_entry["7_mm"] = seven_mm
            dict_entry["8_mm"] = eight_plus_mm
            dict_entry["8_mm_p_fail"] = eight_plus_fail
            dict_entry["failure"] = fail_tally
            dict_entry["pos_elims"] = pos_elims
            dict_entry["neg_elims"] = neg_elims
            dict_entry["recall_TP"] = "match_1_2"
            dict_entry["recall_FN"] = "3_mm_p_fail"

            ''' append dict entry to list '''
            data_list.append(dict_entry)

    ''' Finish up with timestamp '''
    table_dict["Timestamp"] = str(datetime.now())

    return table_dict


def get_arguments():
    parser = ap.ArgumentParser(prog='parse.py',
                               description="""Parse""")
    parser.add_argument('-r', '--resource_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="directory for resource files")
    parser.add_argument('-R', '--results_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="directory for results files")
    parser.add_argument('-f', '--fasta_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="location of genome fasta files")
    return parser.parse_args()



def main():
    ''' Get command line argument '''
    args = get_arguments()

    resource_dir = args.resource_directory[0]
    results_dir = args.results_directory[0]
    fasta_dir = args.fasta_directory[0]

    # File Names
    meta_filename = os.path.join(resource_dir, "nextstrain_ncov_metadata.txt")
    sum_data_file = os.path.join(resource_dir, "esummary_ncbi.txt")
    output_file = os.path.join(results_dir, "summary_table.json")
    short_assay_list_file = os.path.join(resource_dir, "reduced_assays.txt")

    ''' Initialize table dictionary '''
    table_dict = { "data": [ ] }

    ''' Get short assay list '''
    short_assay_list = []
    with open(short_assay_list_file, 'r') as read_file:
        for line in read_file.readlines():
            fields = line.split()
            assay_label = fields[0]
            short_assay_list.append(assay_label)


    ''' Find most recent json file output from Assay Monitor '''
    json_file_list = [os.path.join(results_dir, f) for f in os.listdir(results_dir)
                            if re.match(r'2020_.*Assay_Results\.json$', f) ]
    json_file_list.sort()
    json_file = json_file_list[-1]

    ''' Create Summary table '''
    print()
    print("Creating summary table json from:", json_file)
    table_dict = generate_table(json_file, meta_filename, sum_data_file, table_dict, short_assay_list)

    ''' Output results as json file '''
    with open(output_file, 'w') as file_handle:
        print(json.dumps(table_dict, sort_keys=True, indent=4), file=file_handle)


    print()
    print("Finished")


if __name__ == "__main__":
    main()

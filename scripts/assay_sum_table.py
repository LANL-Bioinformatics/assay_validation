#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:36:54 2020

@author: adanm
"""

import json
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


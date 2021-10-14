#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:36:54 2020

@author: adanm
"""
from datetime import datetime
import multiprocessing as mp

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

def set_queues_for_table_gen(TP_list_generator, FN_list_generator, procnum):
    pool = mp.Pool(processes=procnum)
    pos_list = pool.map_async(positives_stats, TP_list_generator).get()
    neg_list = pool.map_async(negatives_stats, FN_list_generator).get()
    pos_dict = {
        tuple[0]: tuple[1]
        for tuple in pos_list
    }
    neg_dict = {
        tuple[0]: tuple[1]
        for tuple in neg_list
    }
    return pos_dict, neg_dict


def generate_table_parallel(full_dict, procnum):

    print("\nCreating summary table from json source:")

    # loop through all assays in assay results dict, making generator
    TP_list_generator = (  # of tuples with
        (assay_dict["Name"],  # Assay name
        assay_dict["True Positives"])  # list of True Positives
        for assay_dict in full_dict[assay]
    )


    # loop through all assays in assay results dict, making generator
    FN_list_generator = (  # of tuples containing
        (assay_dict["Name"],  # Assay name
        assay_dict["False Negatives"])  # list of False Negatives
        for assay_dict in full_dict[assay]
    )

    # Generate tallies of TPs and FNs from each assay in parallel
    pos_dict, neg_dict = set_queues_for_table_gen(
            TP_list_generator, FN_list_generator, procnum)


    # Populate list with summary data from each assay
    data_list = []
    for asy_name,TP_tally in pos_dict.items():
        stats_dict = calculate_recall(asy_name, TP_tally, neg_dict[asy_name])
        data_list.append(stats_dict)
    
    # Add data to table dictionary
    table_dict = { "data": data_list }

    # Add timestamp
    table_dict["Timestamp"] = str(datetime.now())

    return table_dict


def positives_stats(argument):
    this_assay, TP_list = argument

    # initialize tallies
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

    # Get mismatch tallies from true positive list
    for seq_entry in TP_list:
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
    
    return (this_assay, [succeed_tally, succeed_one_two, one_mm, two_mm, 
                        two_plus, three_mm, three_plus_fail, four_mm, five_mm, 
                        six_mm, seven_mm, eight_plus_mm, eight_plus_fail, 
                        fail_tally, pos_elims, neg_elims])


def negatives_stats(argument):
    this_assay, FN_list = argument

    # initialize tallies
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

    # Get tallies of false negatives
    for _ in FN_list:
        try:
            three_plus_fail += 1
            eight_plus_fail += 1
            fail_tally += 1
        except:
            neg_elims += 1

    return (this_assay, [succeed_tally, succeed_one_two, one_mm, two_mm, 
                        two_plus, three_mm, three_plus_fail, four_mm, five_mm, 
                        six_mm, seven_mm, eight_plus_mm, eight_plus_fail, 
                        fail_tally, pos_elims, neg_elims])


def calculate_recall(this_assay, TP_tallies, FN_tallies):

    # Sum TP tallies with FN tallies
    total_tallies = []
    for i in range(len(TP_tallies)):
        total_tallies.append(TP_tallies[i] + FN_tallies[i])
    
    # Extract all separate totals
    [   succeed_tally, succeed_one_two, one_mm, two_mm, two_plus, three_mm, 
        three_plus_fail, four_mm, five_mm, six_mm, seven_mm, eight_plus_mm, 
        eight_plus_fail, fail_tally, pos_elims, neg_elims ] = total_tallies

    # Calculate recall
    try:
        TPR = str(succeed_one_two/(succeed_one_two + three_plus_fail))
    except:
        print("this_assay", this_assay)
        print("succeed_one_two", succeed_one_two)
        print("three_plus_fail", three_plus_fail)


    # Prepare dictionary entry
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

    return dict_entry


def generate_table(full_dict):

    print("\nCreating summary table from json source:")

    ''' Initialize table dictionary '''
    table_dict = { "data": [ ] }

    ''' Get data list from table dictionary '''
    data_list = table_dict["data"]

    ''' loop through all assays in assay monitor dictionary '''
    for assay_dict in full_dict[assay]:
        this_assay = assay_dict[assay_name]


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


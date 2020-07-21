#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:19:40 2020

@author: adanm
"""

import re
import os
import json
from glob import glob
import argparse as ap
import pandas as pd


# Dictionary Labels
assay = "Assay"
assay_name = "Name"
FP = "Forward Primer"
RP = "Reverse Primer"
Pr = "Probe"
gp = "gaps"
mm = "mismatches"
positives_type = "True Positives"
negatives_type = "False Negatives"


def parse_tnt_results(full_dictionary, tnt_result):
    tnt_result_filename =tnt_result.split('/')[-1]
    tnt_result_prefix = tnt_result_filename.split('.')[0]
    seq_accession = re.sub(r'_results', '', tnt_result_prefix)

    ''' Parse tnt file for positive assays and relevant info '''
    assay_list = full_dictionary[assay]
    with open(tnt_result, 'r') as file_handle:
        for line in file_handle.readlines():
            fields = line.strip().split(" ")
            if len(fields) >= 3:

                if fields[0] == "name":
                    name = fields[2]
#                    print("#### name = {}".format(name))
                    if not assay_list:
#                        print("assay_list empty: Creating assay entry in list")
                        assay_dict = { assay_name : name }
                        assay_dict.update({positives_type : [] })
                        assay_dict.update({negatives_type : [] })
                        assay_list.append(assay_dict)

                    else:
#                        print("assay_list not empty: Find assay in list")
                        if name in [ dictionary[assay_name] for dictionary in assay_list ]:
                            for dictionary in assay_list:
                                if dictionary[assay_name] == name:
                                    assay_dict = dictionary
                        else:
#                            print("assay not in list: Creating assay entry in lis")
                            assay_dict = { assay_name : name }
                            assay_dict.update({positives_type : [] })
                            assay_dict.update({negatives_type : [] })
                            assay_list.append(assay_dict)

#                    print("First time though: Create sequence accesion entry for", seq_accession)

                    sequence_dict = { "Accession" : seq_accession, "Values" : {},
                                      "Alignments" : {}, "Sequences" : {} ,
                                      "Composition" : {}, "Heuristics" : {},
                                      "Thermo" : {} }
                    positives_list = assay_dict[positives_type]
                    sequence_dict["Values"].update({"Forward Primer" : {} })
                    sequence_dict["Values"].update({"Reverse Primer" : {} })
                    sequence_dict["Values"].update({"Probe" : {} })
                    sequence_dict["Alignments"].update({FP: {}})
                    sequence_dict["Alignments"].update({RP: {}})
                    sequence_dict["Alignments"].update({Pr: {}})

                ''' Populate mismatches and gaps fields '''
                if fields[2] == "mismatches":
                    if fields[0] == "reverse":
                        sequence_dict["Values"]["Reverse Primer"].update({"mismatches" : fields[4]})
                    if fields[0] == "forward":
                        sequence_dict["Values"]["Forward Primer"].update({"mismatches" : fields[4]})
                if fields[1] == "mismatches":
                    sequence_dict["Values"]["Probe"].update({"mismatches" : fields[3]})
                if fields[2] == "gaps":
                    if fields[0] == "reverse":
                        sequence_dict["Values"]["Reverse Primer"].update({"gaps" : fields[4]})
                    if fields[0] == "forward":
                        sequence_dict["Values"]["Forward Primer"].update({"gaps" : fields[4]})
                if fields[1] == "gaps":
                    sequence_dict["Values"]["Probe"].update({"gaps" : fields[3]})


                ''' Populate Alignments fields'''
                if fields[2] == "align":
                    if fields[0] == "forward" and fields[1] == "primer":
                        if fields[3] == "":
                            sequence_dict["Alignments"][FP].update({"pairing": " ".join(fields[6:])})
                        elif fields[3] == "dimer":
                            sequence_dict["Alignments"][FP].update({"dimer alignment size": fields[7]})
                        else:
                            sequence_dict["Alignments"][FP].update({fields[3].split('-')[0]: fields[3]})


                    if fields[0] == "reverse" and fields[1] == "primer":
                        if fields[3] == "":
                            sequence_dict["Alignments"][RP].update({"pairing": " ".join(fields[6:])})
                        elif fields[3] == "dimer":
                            sequence_dict["Alignments"][RP].update({"dimer alignment size": fields[7]})
                        else:
                            sequence_dict["Alignments"][RP].update({fields[3].split('-')[0]: fields[3]})

                if fields[0] == "probe"  and fields[1] == "align":
                    if fields[2] == "":
                        sequence_dict["Alignments"][Pr].update({"pairing": " ".join(fields[5:])})
                    elif fields[2] == "dimer":
                        sequence_dict["Alignments"][Pr].update({"dimer alignment size": fields[6]})
                    else:
                        sequence_dict["Alignments"][Pr].update({fields[2].split('-')[0]: fields[2]})


                ''' Populate Sequences fields'''
                if fields[1] == "primer" and fields[2] == "=":
                    if fields[0] == "forward":
                        sequence_dict["Sequences"].update({FP: fields[3]})
                    if fields[0]== "reverse":
                        sequence_dict["Sequences"].update({RP: fields[3]})
                if fields[0] == "probe" and fields[1] == "=":
                    sequence_dict["Sequences"].update({Pr: fields[2]})


                ''' Populate melting temp fields '''
                if fields[2] == "tm":
                    if fields[0] == "reverse":
                        sequence_dict["Values"][RP].update({"tm" : fields[4]})
                    if fields[0] == "forward":
                        sequence_dict["Values"][FP].update({"tm" : fields[4]})
                if fields[2] == "hairpin" and fields[3] == "tm":
                    if fields[0] == "reverse":
                        sequence_dict["Values"][RP].update({"hairpin tm" : fields[5]})
                    if fields[0] == "forward":
                        sequence_dict["Values"][FP].update({"hairpin tm" : fields[5]})
                if fields[2] == "homodimer" and fields[3] == "tm":
                    if fields[0] == "reverse":
                        sequence_dict["Values"][RP].update({"homodimer tm" : fields[5]})
                    if fields[0] == "forward":
                        sequence_dict["Values"][FP].update({"homodimer tm" : fields[5]})
                if fields[1] == "tm":
                    sequence_dict["Values"][Pr].update({"tm" : fields[3]})
                if fields[1] == "hairpin" and fields[2] == "tm":
                    sequence_dict["Values"][Pr].update({"hairpin tm" : fields[4]})
                if fields[1] == "homodimer" and fields[2] == "tm":
                    sequence_dict["Values"][Pr].update({"homodimer tm" : fields[4]})
                if fields[0] == "heterodimer" and fields[1] == "tm":
                    sequence_dict["Values"].update( {"Heterodimer" : {"heterodimer tm" : fields[3]} } )


                ''' Populate Composition fields '''
                if fields[1] == "3'" and fields[2] == "clamp":
                    if fields[0] == "min":
                        sequence_dict["Composition"].update( { "min 3' clamp" : fields[4]})
                    if fields[0] == "max":
                        sequence_dict["Composition"].update( { "max 3' clamp" : fields[4]})
                if fields[1] == "primer" and fields[2] == "%GC":
                    if fields[0] == "forward":
                        sequence_dict["Composition"].update( { "forward primer %GC " : fields[4]} )
                    if fields[0] == "reverse":
                        sequence_dict["Composition"].update( { "reverse primer %GC " : fields[4]} )
                if fields[0] == "amplicon" and fields[1] == "range":
                    sequence_dict["Composition"].update( { "amplicon range" : [ fields[3], fields[5]]})
                if fields[0] == "amplicon" and fields[1] == "length":
                    sequence_dict["Composition"].update( { "amplicon length" : fields[3]})
                if fields[0] == "Forward" and fields[1] == "primer":
                    sequence_dict["Composition"].update( { "primer contained message" : line.strip()})
                if fields[0] == "probe" and fields[1] == "%GC":
                    sequence_dict["Composition"].update( { "probe %GC" : fields[3]} )
                if fields[0] == "probe" and fields[1] == "range":
                    sequence_dict["Composition"].update( { "probe range" : [ fields[3], fields[5]] })
                if fields[0] == "probe" and fields[1] == "contained":
                    sequence_dict["Composition"].update( { "probe contained message" : line.strip()})


                ''' Populate Heuristics fields '''
                if fields[2] == "heuristics":
                    if fields[0] == "forward":
                        sequence_dict["Heuristics"].update( { "forward primer heuristics" : fields[4] } )
                    if fields[0] == "reverse":
                        sequence_dict["Heuristics"].update( { "reverse primer heuristics" : fields[4] } )


                ''' Populate Thermo Fields '''
                if fields[2].split("[")[0] == "dG":
                    if fields[0] == "forward":
                        sequence_dict["Thermo"].update( { FP : {} } )
                        sequence_dict["Thermo"][FP].update( { "dG" : fields[2].split("[")[1][0:-1] } )
                        sequence_dict["Thermo"][FP].update( { "dH" : fields[4].split("[")[1][0:-1] } )
                        sequence_dict["Thermo"][FP].update( { "dS" : fields[6].split("[")[1][0:-1] } )
                    if fields[0] == "reverse":
                        sequence_dict["Thermo"].update( { RP : {} } )
                        sequence_dict["Thermo"][RP].update( { "dG" : fields[2].split("[")[1][0:-1] } )
                        sequence_dict["Thermo"][RP].update( { "dH" : fields[4].split("[")[1][0:-1] } )
                        sequence_dict["Thermo"][RP].update( { "dS" : fields[6].split("[")[1][0:-1] } )
                if fields[1].split("[")[0] == "dG":
                    sequence_dict["Thermo"].update( { Pr : {} } )
                    sequence_dict["Thermo"][Pr].update( { "dG" : fields[1].split("[")[1][0:-1] } )
                    sequence_dict["Thermo"][Pr].update( { "dH" : fields[3].split("[")[1][0:-1] } )
                    sequence_dict["Thermo"][Pr].update( { "dS" : fields[5].split("[")[1][0:-1] } )


            ''' Add common name for genome '''
            if line[0] == ">":
                sequence_dict["Common Name"] = line[1:].strip().split('|')[0]
                sequence_dict["Accession"] = line.split('|')[1]

                positives_list.append(sequence_dict)

    return full_dictionary



def determine_negatives(assay_file, tnt_result, full_dictionary, positives_type, negatives_type):
    tnt_result_filename =tnt_result.split('/')[-1]
    tnt_result_prefix = tnt_result_filename.split('.')[0]
    seq_accession = re.sub(r'_results', '', tnt_result_prefix)

    ''' Negatives '''
    # Find which assays were negative from original list
    assay_list = full_dictionary[assay]
    for assay_dictionary in assay_list:
#        name = assay_dictionary["Name"]
        positives_list = assay_dictionary[positives_type]
        negatives_list = assay_dictionary[negatives_type]
        if not seq_accession in [ acc_dictionary['Accession'] for acc_dictionary in positives_list]:
            sequence_dict = { "Accession" : seq_accession}
            negatives_list.append(sequence_dict)
    return full_dictionary


def get_working_assays(assay_dict, positives_type, negatives_type):
    positives = []
    negatives = []
    pos_acc_list = assay_dict[positives_type]
    neg_acc_list = assay_dict[negatives_type]
    for acc_entry in pos_acc_list:
        accession = acc_entry["Accession"]
        positives.append(acc_entry["Accession"])
    for acc_entry in neg_acc_list:
        accession = acc_entry["Accession"]
        negatives.append(acc_entry["Accession"])
    return positives, negatives


def create_match_tables_reduced(in_file, results_path, seq_file, short_assay_list_file):

    ''' Get short assay list '''
    short_assay_list = []
    with open(short_assay_list_file, 'r') as read_file:
        for line in read_file.readlines():
            fields = line.split()
            assay_label = fields[0]
            short_assay_list.append(assay_label)


    seq_list = []
    with open(seq_file, 'r') as file_handle:
        lines = file_handle.readlines()
        for line in lines:
            seq_list.append(line.strip())

    length = len(seq_list)
    dict_of_success = {}
    dict_of_fail = {}
    dict_of_both = {}
#    json_file_path = os.path.join(results_path, in_file)
    with open(in_file) as json_file:
        full_dict = json.load(json_file)
        for assay_dict in full_dict[assay]:
            this_assay = assay_dict[assay_name]
            if this_assay in short_assay_list:

                positives, negatives = get_working_assays(assay_dict, "True Positives", "False Negatives")
                table_line = [None]*(length)
                for acc in positives:
                    try:
                        idx = seq_list.index(acc)
                        mismatches = 0
                        for seq_entry in assay_dict["True Positives"]:
                            if seq_entry["Accession"] == acc:

                                mismatches_fp = int(seq_entry["Values"]["Forward Primer"]['mismatches'])
                                if 'mismatches' in seq_entry["Values"]["Probe"].keys():
                                    mismatches_pr = int(seq_entry["Values"]["Probe"]['mismatches'])
                                else:
                                    mismatches_pr = 0
                                mismatches_rp = int(seq_entry["Values"]["Reverse Primer"]['mismatches'])
                                mismatches = max(mismatches_fp, mismatches_pr, mismatches_rp)

    #                    print(idx, seq_list[idx])
                        table_line[idx] = mismatches
                    except:
                        with open("tnt_pass_tally_file.txt", 'a') as err_file:
                            print("pass pos", file=err_file)


                for acc in negatives:
                    try:
                        idx = seq_list.index(acc)
    #                    print(idx, seq_list[idx])
                        table_line[idx] = -1
                    except:
                        with open("tnt_pass_tally_file.txt", 'a') as err_file:
                            print("pass neg", file=err_file)
    #            print(table_line)
                dict_of_success.update({this_assay : positives})
                dict_of_fail.update({this_assay : negatives})
                dict_of_both.update({this_assay : table_line})

    with open(os.path.join(results_path, "match_table.csv"), 'w') as file_handle:
        text_line = " ," + ",".join([ str(item) for item in seq_list ])
        print(text_line, file=file_handle)

        for key,val in dict_of_both.items():

            text_line = key + "," + ",".join([ str(item) for item in val ])
            print(text_line, file=file_handle)

    print("Wrote file: {}match_table.csv".format(results_path))


def make_sequence_list_file(seq_file, target_fna_file_path):
    fna_filepath_list = glob(os.path.join(target_fna_file_path, "*.fna"))
    fna_file_list = [ fna_file.split('/')[-1] for fna_file in fna_filepath_list ]
    accession_list = [ each_file.split('.')[0] for each_file in fna_file_list ]
    accession_list.sort()
    with open(seq_file, 'w') as write_handle:
        for each in accession_list:
            print(each, file=write_handle)


def give_complement_base(base):
    letter = base.upper()
    if letter=='A':
        return 'T'
    elif letter=='T':
        return 'A'
    elif letter=='G':
        return 'C'
    elif letter=='C':
        return 'G'
    else:
        return "N"


def reverse_complement(sequence):
    return ''.join([give_complement_base(base) for base in reversed(sequence)])


def complement(sequence):
    return ''.join([give_complement_base(base) for base in sequence])


def get_del_ct(primer_bases, target_bases, three_prime_table):
    table = pd.read_table(three_prime_table, index_col='Row')
    return table.loc[primer_bases, target_bases]


def is_fatal_terminal_mismatch(primer_bases, target_bases, three_prime_table, del_ct_thresh):
    primer_bases, target_bases = process_Ns(primer_bases, target_bases)
    result = get_del_ct(primer_bases, target_bases, three_prime_table)
    if result != 0.0:
        print(result)
    return get_del_ct(primer_bases, target_bases, three_prime_table) > del_ct_thresh


def process_Ns(primer_bases, target_bases):
    if "N" in primer_bases:
        print(primer_bases, target_bases)
        if primer_bases[0] == "N":
            primer_bases = "{}{}".format(target_bases[0], primer_bases[1])
        if primer_bases[1] == "N":
            primer_bases = "{}{}".format(primer_bases[0], target_bases[1])
        print(primer_bases, target_bases)

    if "N" in target_bases:
        print(primer_bases, target_bases)
        if target_bases[0] == "N":
            target_bases = "{}{}".format(primer_bases[0], target_bases[1])
        if target_bases[1] == "N":
            target_bases = "{}{}".format(target_bases[0], primer_bases[1])
        print(primer_bases, target_bases)
    return primer_bases, target_bases


def filter_three_prime(Full_Dict, three_prime_table, del_ct_thresh):
    assay_dict_list = Full_Dict["Assay"]
    for each_assay in assay_dict_list:
        FN_list = each_assay["False Negatives"]
        TP_list = each_assay["True Positives"]
        new_TP_list = []
        for each_pos in TP_list:

            ''' Get last two bases of Primers '''
            # Forward Primer
            for_prim = each_pos["Alignments"]["Forward Primer"]["5'"]
            for_seq = for_prim.split('-')[1][-2:]

            # Reverse Primer
            rev_prim = each_pos["Alignments"]["Reverse Primer"]["5'"]
            rev_seq = rev_prim.split('-')[1][-2:]

            ''' Get corresponding bases from Targets '''
            # Forward Primer
            for_prim_targ = each_pos["Alignments"]["Forward Primer"]["3'"]
            for_targ_seq = for_prim_targ.split('-')[1][-2:]
            for_targ_seq = complement(for_targ_seq)

            # Reverse Primer
            rev_prim_targ = each_pos["Alignments"]["Reverse Primer"]["3'"]
            rev_targ_seq = rev_prim_targ.split('-')[1][-2:]
            rev_targ_seq = complement(rev_targ_seq)

            ''' If terminal MM put in False Negative List; else put in new TP list '''
            if is_fatal_terminal_mismatch(for_seq, for_targ_seq, three_prime_table, del_ct_thresh) or \
                is_fatal_terminal_mismatch(rev_seq, rev_targ_seq, three_prime_table, del_ct_thresh):
                    print("fatal")
                    print(each_assay["Name"])
                    print(each_pos["Accession"])
                    print(for_seq, for_targ_seq, rev_seq, rev_targ_seq)
                    print()
                    FN_list.append(each_pos)
            else:
                new_TP_list.append(each_pos)

        each_assay["True Positives"] = new_TP_list

    return Full_Dict



def get_arguments():
    parser = ap.ArgumentParser(prog='assay_monitor.py', description="""Assay Monitor""")
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
    ''' Get command line arguments '''
    args = get_arguments()

    resource_dir = args.resource_directory[0]
    results_dir = args.results_directory[0]
    fasta_dir = args.fasta_directory[0]
    tnt_res_dir = args.tnt_results_directory[0]

    # File Names
    output_file_path = os.path.join(results_dir, "Assay_Results.json")
    assay_file = os.path.join(resource_dir, "assays.txt")

    # Initialize assay dictionary
    Full_Dict = {assay : []}

    ''' Parse tnt output for True Positives into multi-level dictionary for json file: '''
    tnt_results_list = [ os.path.join(tnt_res_dir, f) for f in os.listdir(tnt_res_dir) if re.match(r'.*_results\.out', f) ]
    for tnt_result in tnt_results_list:
        Full_Dict = parse_tnt_results(Full_Dict, tnt_result)


    ''' False Negatives'''
    ''' Parse and Print false negative results '''
    for tnt_result in tnt_results_list:
        Full_Dict = determine_negatives(assay_file, tnt_result, Full_Dict, "True Positives", "False Negatives")


    with open(output_file_path, 'w') as file_handle:
        print(json.dumps(Full_Dict, sort_keys=True, indent=4), file=file_handle)


    print()
    print("Finished")


if __name__ == "__main__":
    main()

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
import pandas as pd
from datetime import datetime
import multiprocessing as mp


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


def get_isolate(proto_isolate):
    isolate = proto_isolate.strip().strip('>')
    return ''.join(isolate.split(' '))



def get_acc_list(file_ns):
    metafile = file_ns.metafile
    use_gisaid = file_ns.use_gisaid
    accessions_file = file_ns.accession_file

    iso_acc_dict = {}

    if use_gisaid:
        with open(metafile, 'r') as read_meta:
            # Get column headers for "field names"
            column_headers = read_meta.readline().strip().split('\t')

            # Read through remainder of metafile
            for metaline in read_meta.readlines():
                # Assign each column value to appropriate field name in dict
                metaline_fields = metaline.strip().split('\t')
                meta_dict = {
                    column_headers[i]: metaline_fields[i]
                    for i in range(len(metaline_fields))
                }

                # Add iso-acc pair to dict
                isolate = meta_dict["Virus name"]
                accession = meta_dict["Accession ID"]
                iso_acc_dict[isolate] = accession

                # Add iso-acc pair for isolates with no spaces
                spaceless_isolate = re.sub(r"\s+", "", isolate)
                iso_acc_dict[spaceless_isolate] = accession
    else:
        with open(accessions_file, 'r') as read_accs:
            for acc_line in read_accs.readlines():
                accession = acc_line.strip()
                iso_acc_dict[accession] = accession
   
    return iso_acc_dict
            
            
    
    
def return_accession(proto_isolate, iso_acc_dict, iso_idx):

    pre_isolate = get_isolate(proto_isolate)
    isolate_fields = pre_isolate.strip().split('|')
    isolate = isolate_fields[iso_idx]
    accession = iso_acc_dict[isolate]

    try:

        return accession

    except:
        print(isolate, "not found in metafile")
        return None
    
    
    
    
def get_header(isolate, metafile):

    header_dict = {}

    with open(metafile, 'r') as read_meta:
        column_headers = read_meta.readline().strip().split('\t')

        for metaline in read_meta.readlines():
            metaline_fields = metaline.strip().split('\t')
            # print(metaline_fields)
            if isolate in metaline_fields:
                for i in range(len(column_headers)):
                    header_dict[column_headers[i]] = metaline_fields[i]

                try:

                    isolate = header_dict["strain"]
                    virus = header_dict["virus"]
                    accession = header_dict["gisaid_epi_isl"]
                    collection_date = header_dict["date"]
                    region = header_dict["region"]
                    country = header_dict["country"]
                    division = header_dict["division"]
                    location = header_dict["location"]
                    region_exposure = header_dict["region_exposure"]
                    country_exposure = header_dict["country_exposure"]
                    division_exposure = header_dict["division_exposure"]
                    segment = header_dict["segment"]
                    length = header_dict["length"]
                    host = header_dict["host"]
                    age = header_dict["age"]
                    sex = header_dict["sex"]
                    originating_lab = header_dict["originating_lab"]
                    submitting_lab = header_dict["submitting_lab"]
                    authors = header_dict["authors"]
                    url = header_dict["url"]
                    title  = header_dict["title"]
                    paper_url = header_dict["paper_url"]
                    date_submitted = header_dict["date_submitted"]
                    break
                except:
                    print(isolate, "######################################################")

    try:
        if virus == "ncov":
            viral_prefix = "hCoV-19/"
        else:
            viral_prefix = ""

        header = ">{}{}|{}|{}|{}|{}: {}|{}|{}|{}|{}/{}/{}|{}/{}/{}/{}|{}|{}|{}|{}|{}|{}|{}|{}".format(viral_prefix, isolate, accession, collection_date, date_submitted, country, division,length, "", "", region, country, division, location, region_exposure, country_exposure, division_exposure, segment, host, age, sex, originating_lab, submitting_lab, authors, title)
        # header = ">{}{}|{}|{}".format(viral_prefix, isolate, accession, collection_date)

        return accession, header

    except:
        print(isolate, "not found in metafile")
        return None, None
    
    
    
def split_gisaid_file_with_meta(gisaid_file, download_dir, metafile, limit_gisaid):
    ''' ################################################ XXXX ### '''

    with open(gisaid_file, 'r') as read_gisaid:
        line = read_gisaid.readline()
        isolate = get_isolate(line)
        line2 = []
        accession, header = get_header(isolate, metafile)

        if not accession == None:
            filename = accession + ".fna"
            output_path = os.path.join(download_dir, filename)
            with open(output_path, "w") as write_handle:
                print(header, file=write_handle)
            # print(line.strip())


        for line in read_gisaid.readlines():
            if line[0] == ">":
                # print(line.strip())
                line2_str = "".join(line2)

                if not accession == None:
                    with open(output_path, "a") as write_handle:
                        print(line2_str.strip(), file=write_handle)

                isolate = get_isolate(line)
                line2 = []
                accession, header = get_header(isolate, metafile)

                if not accession == None:
                    filename = accession + ".fna"
                    output_path = os.path.join(download_dir, filename)
                    with open(output_path, "w") as write_handle:
                        print(header, file=write_handle)
            else:
                # print(line[0:40])
                line2.append(line.strip())


        line2_str = "".join(line2)

        if not accession == None:
            with open(output_path, "a") as write_handle:
                print(line2_str.strip(), file=write_handle)


def next_fields(tnt_lines):
    line = next(tnt_lines)
    fields = line.strip().split(" ")
    return fields


def fill_seq_dict(tnt_lines, iso_acc_dict, iso_idx):

    # Initialize Sequence Dict
    sequence_dict = { "Values" : {},
                        "Alignments" : {}, "Sequences" : {} ,
                        "Composition" : {}, "Heuristics" : {},
                        "Thermo" : {} }
    sequence_dict["Values"].update({"Forward Primer" : {} })
    sequence_dict["Values"].update({"Reverse Primer" : {} })
    sequence_dict["Values"].update({"Probe" : {} })
    sequence_dict["Alignments"].update({FP: {}})
    sequence_dict["Alignments"].update({RP: {}})
    sequence_dict["Alignments"].update({Pr: {}})

    ''' Populate Sequences fields'''
    fields = next_fields(tnt_lines)
    sequence_dict["Sequences"].update({FP: fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Sequences"].update({RP: fields[4]})

    ''' Populate melting temp fields '''
    fields = next_fields(tnt_lines)
    sequence_dict["Values"][FP].update({"tm" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"][RP].update({"tm" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"][FP].update({"hairpin tm" : fields[5]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"][RP].update({"hairpin tm" : fields[5]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"][FP].update({"homodimer tm" : fields[5]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"][RP].update({"homodimer tm" : fields[5]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"].update( {"Heterodimer" : {"heterodimer tm" : fields[3]} } )

    ''' Populate Thermo Fields '''
    fields = next_fields(tnt_lines)
    sequence_dict["Thermo"].update( { FP : {} } )
    sequence_dict["Thermo"][FP].update( { "dG" : fields[2].split("[")[1][0:-1] } )
    sequence_dict["Thermo"][FP].update( { "dH" : fields[4].split("[")[1][0:-1] } )
    sequence_dict["Thermo"][FP].update( { "dS" : fields[6].split("[")[1][0:-1] } )
    fields = next_fields(tnt_lines)
    sequence_dict["Thermo"].update( { RP : {} } )
    sequence_dict["Thermo"][RP].update( { "dG" : fields[2].split("[")[1][0:-1] } )
    sequence_dict["Thermo"][RP].update( { "dH" : fields[4].split("[")[1][0:-1] } )
    sequence_dict["Thermo"][RP].update( { "dS" : fields[6].split("[")[1][0:-1] } )

    ''' Populate mismatches and gaps fields '''
    fields = next_fields(tnt_lines)
    sequence_dict["Values"]["Forward Primer"].update({"mismatches" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"]["Reverse Primer"].update({"mismatches" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"]["Forward Primer"].update({"gaps" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Values"]["Reverse Primer"].update({"gaps" : fields[4]})

    ''' Populate Composition fields '''
    fields = next_fields(tnt_lines)
    sequence_dict["Composition"].update( { "min 3' clamp" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Composition"].update( { "max 3' clamp" : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Composition"].update( { "forward primer %GC " : fields[4]} )
    fields = next_fields(tnt_lines)
    sequence_dict["Composition"].update( { "reverse primer %GC " : fields[4]} )

    ''' Populate Heuristics fields '''
    fields = next_fields(tnt_lines)
    sequence_dict["Heuristics"].update( { "forward primer heuristics" : fields[4] } )
    fields = next_fields(tnt_lines)
    sequence_dict["Heuristics"].update( { "reverse primer heuristics" : fields[4] } )

    ''' Populate Composition fields '''
    fields = next_fields(tnt_lines)
    sequence_dict["Composition"].update( { "amplicon range" : [ fields[3], fields[5]]})
    fields = next_fields(tnt_lines)
    sequence_dict["Composition"].update( { "amplicon length" : fields[3]})
    line = next(tnt_lines)
    sequence_dict["Composition"].update( { "primer contained message" : line.strip()})

    ''' Probe Fields '''
    fields = next_fields(tnt_lines)
    if fields[0] == "probe":
        sequence_dict["Sequences"].update({Pr: fields[3]})
        fields = next_fields(tnt_lines)
        sequence_dict["Values"][Pr].update({"tm" : fields[3]})
        fields = next_fields(tnt_lines)
        sequence_dict["Values"][Pr].update({"hairpin tm" : fields[4]})
        fields = next_fields(tnt_lines)
        sequence_dict["Values"][Pr].update({"homodimer tm" : fields[4]})
        fields = next_fields(tnt_lines)
        sequence_dict["Thermo"].update( { Pr : {} } )
        sequence_dict["Thermo"][Pr].update( { "dG" : fields[1].split("[")[1][0:-1] } )
        sequence_dict["Thermo"][Pr].update( { "dH" : fields[3].split("[")[1][0:-1] } )
        sequence_dict["Thermo"][Pr].update( { "dS" : fields[5].split("[")[1][0:-1] } )
        fields = next_fields(tnt_lines)
        sequence_dict["Values"]["Probe"].update({"mismatches" : fields[3]})
        fields = next_fields(tnt_lines)
        sequence_dict["Values"]["Probe"].update({"gaps" : fields[3]})
        fields = next_fields(tnt_lines)
        sequence_dict["Composition"].update( { "probe %GC" : fields[3]} )
        fields = next_fields(tnt_lines)
        sequence_dict["Composition"].update( { "probe range" : [ fields[3], fields[5]] })
        line = next(tnt_lines)
        sequence_dict["Composition"].update( { "probe contained message" : line.strip()})
        fields = next_fields(tnt_lines)

    ''' Populate Alignments fields'''
    sequence_dict["Alignments"][FP].update({fields[3] : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][FP].update({"pairing": " ".join(fields[6:])})
    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][FP].update({fields[3] : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][FP].update({"dimer alignment size": fields[7]})

    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][RP].update({fields[3] : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][RP].update({"pairing": " ".join(fields[6:])})
    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][RP].update({fields[3] : fields[4]})
    fields = next_fields(tnt_lines)
    sequence_dict["Alignments"][RP].update({"dimer alignment size": fields[7]})

    line = next(tnt_lines)
    fields = line.strip().split(" ")
    if fields[0] == "probe":
        sequence_dict["Alignments"][Pr].update({fields[2] : fields[3]})
        fields = next_fields(tnt_lines)
        sequence_dict["Alignments"][Pr].update({"pairing": " ".join(fields[5:])})
        fields = next_fields(tnt_lines)
        sequence_dict["Alignments"][Pr].update({fields[2] : fields[3]})
        fields = next_fields(tnt_lines)
        sequence_dict["Alignments"][Pr].update({"dimer alignment size": fields[6]})
        line = next(tnt_lines)

    ''' Add common name for genome '''
    sequence_dict["Common Name"] = line[1:].strip()
    accession = return_accession(
        line[1:].strip(), iso_acc_dict, iso_idx)
    sequence_dict["Accession"] = accession

    ''' Run out remainder of lines for this entry '''
    line = next(tnt_lines)
    line = next(tnt_lines)

    return sequence_dict


def assay_pos_list(tnt_lines, iso_acc_dict, iso_idx):
    positives_list = []
    while line := next(tnt_lines, False):
        if line.strip() == ("###########################################"
                            + "##########################################"):
            break
        fields = line.strip().split(" ")
        name = fields[2]
        
        sequence_dict = fill_seq_dict(tnt_lines, iso_acc_dict, iso_idx)

        ''' Add this sequence dict to positives list '''
        positives_list.append(sequence_dict)
    return positives_list, name, line


def make_assay_list(tnt_lines, iso_acc_dict, iso_idx):

    # Get first line and discard
    line = next(tnt_lines)
    if line.strip() == (
            "###########################################" +
            "##########################################"):
        pass

    # Initialize assay list (for assay dicts)
    assay_list = []

    # Loop while lines exist and create assay dicts, add these to list
    while line:
        positives_list, name, line = assay_pos_list(tnt_lines, iso_acc_dict, iso_idx)

        # Add positives to new assay dict (append to assay list)
        assay_dict = { assay_name : name }
        assay_dict.update({positives_type : positives_list })
        assay_dict.update({negatives_type : [] })
        assay_list.append(assay_dict)

    return assay_list


def brute_parse_tnt_results(tnt_result, iso_acc_dict, iso_idx):
    ''' Parse tnt file for positive assays and relevant info '''

    # Read tnt results file
    with open(tnt_result, 'r') as file_handle:
        tnt_lines_list = file_handle.readlines()

    # Convert tnt file lines to generator
    tnt_lines = (line for line in tnt_lines_list)

    # Create list of assay dicts containing tnt results
    assay_list = make_assay_list(tnt_lines, iso_acc_dict, iso_idx)

    # Add assay list as dictionary entry
    full_dictionary = {assay : assay_list}

    # Add timestamp
    full_dictionary.update({ "Timestamp" : str(datetime.now()) })

    return full_dictionary


def parse_tnt_results(tnt_result, iso_acc_dict, iso_idx):
    ''' Parse tnt file for positive assays and relevant info '''

    # Initialize multilevel dictionary for tnt assay results
    full_dictionary = {assay : []}

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

                    sequence_dict = { "Values" : {},
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
                            sequence_dict["Alignments"][FP].update({fields[3] : fields[4]})


                    if fields[0] == "reverse" and fields[1] == "primer":
                        if fields[3] == "":
                            sequence_dict["Alignments"][RP].update({"pairing": " ".join(fields[6:])})
                        elif fields[3] == "dimer":
                            sequence_dict["Alignments"][RP].update({"dimer alignment size": fields[7]})
                        else:
                            sequence_dict["Alignments"][RP].update({fields[3] : fields[4]})

                if fields[0] == "probe"  and fields[1] == "align":
                    if fields[2] == "":
                        sequence_dict["Alignments"][Pr].update({"pairing": " ".join(fields[5:])})
                    elif fields[2] == "dimer":
                        sequence_dict["Alignments"][Pr].update({"dimer alignment size": fields[6]})
                    else:
                        sequence_dict["Alignments"][Pr].update({fields[2] : fields[3]})

#######################################################
                ''' Populate Sequences fields'''
                if fields[1] == "primer" and fields[2] == "=":
                    if fields[0] == "forward":
                        sequence_dict["Sequences"].update({FP: fields[4]})
                    if fields[0]== "reverse":
                        sequence_dict["Sequences"].update({RP: fields[4]})
                if fields[0] == "probe" and fields[1] == "=":
                    sequence_dict["Sequences"].update({Pr: fields[3]})


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
                sequence_dict["Common Name"] = line[1:].strip()
                accession = return_accession(
                    line[1:].strip(), iso_acc_dict, iso_idx)
                sequence_dict["Accession"] = accession

                positives_list.append(sequence_dict)

    ''' Get current time for time stamp and filename '''
    current_time = datetime.now()

    ''' Add timestamp '''
    current_time_string = str(current_time)
    full_dictionary.update({ "Timestamp" : current_time_string })

    return full_dictionary



def determine_negatives(
        accession_list, full_dictionary, positives_type, negatives_type):

    print("\nDetermining Negatives and Adding To Metadata Dictionary...")

    assay_list = full_dictionary[assay]

    for assay_dictionary in assay_list:
        positives_list = assay_dictionary[positives_type]

        pos_accs = [ acc_dictionary['Accession'] 
            for acc_dictionary in positives_list ]

        neg_accs = set(accession_list) - set(pos_accs)

        neg_dict = [
            {"Accession": neg_acc}
            for neg_acc in neg_accs
        ]
        assay_dictionary[negatives_type] = neg_dict

    # Add timestamp
    full_dictionary.update({ "Timestamp" : str(datetime.now()) })

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
    # if result != 0.0:
    #     print(result)
    return get_del_ct(primer_bases, target_bases, three_prime_table) > del_ct_thresh


def process_Ns(primer_bases, target_bases):
    if "N" in primer_bases:
        # print(primer_bases, target_bases)
        if primer_bases[0] == "N":
            primer_bases = "{}{}".format(target_bases[0], primer_bases[1])
        if primer_bases[1] == "N":
            primer_bases = "{}{}".format(primer_bases[0], target_bases[1])
        # print(primer_bases, target_bases)

    if "N" in target_bases:
        # print(primer_bases, target_bases)
        if target_bases[0] == "N":
            target_bases = "{}{}".format(primer_bases[0], target_bases[1])
        if target_bases[1] == "N":
            target_bases = "{}{}".format(target_bases[0], primer_bases[1])
        # print(primer_bases, target_bases)
    return primer_bases, target_bases


def set_queues_for_three_prime(
        Full_Dict, three_prime_table, del_ct_thresh, procnum):

    print("\nRemoving 3' mismatches over threshold from True Positives...")
    assay_dict_list = Full_Dict["Assay"]

    gen_of_TP_lists = (
        (each_assay["Name"], # Tuple with assay name, and...
            each_assay["True Positives"], # TP_list in each_assay, and...
            three_prime_table, del_ct_thresh # other needed args
        )
        for each_assay in Full_Dict["Assay"] # for each_assay in assay list
    )

    pool = mp.Pool(processes=procnum)
    asy_tuples = pool.map_async(filter_three_prime_parallel, 
                                gen_of_TP_lists).get()
    
    for assay_name,new_TP_list,new_FN_list in asy_tuples:
        for asy in assay_dict_list:
            if asy["Name"] == assay_name:
                asy["True Positives"] = new_TP_list
                asy["False Negatives"].extend(new_FN_list)
                continue

    # Add timestamp
    Full_Dict.update({ "Timestamp" : str(datetime.now()) })

    return Full_Dict


def filter_three_prime_parallel(argument):
    
    assay_name, TP_list, three_prime_table, del_ct_thresh = argument

    new_TP_list = []
    new_FN_list = []
    for each_pos in TP_list:

        ''' Get last two bases of Primers '''
        # Forward Primer
        for_prim = each_pos["Alignments"]["Forward Primer"]["5'"]
        for_seq = for_prim[-2:]

        # Reverse Primer
        rev_prim = each_pos["Alignments"]["Reverse Primer"]["5'"]
        rev_seq = rev_prim[-2:]

        ''' Get corresponding bases from Targets '''
        # Forward Primer
        for_prim_targ = each_pos["Alignments"]["Forward Primer"]["3'"]
        for_targ_seq = for_prim_targ[-2:]
        for_targ_seq = complement(for_targ_seq)

        # Reverse Primer
        rev_prim_targ = each_pos["Alignments"]["Reverse Primer"]["3'"]
        rev_targ_seq = rev_prim_targ[-2:]
        rev_targ_seq = complement(rev_targ_seq)

        ''' If terminal MM put in False Negative List; else put in new TP list '''
        if is_fatal_terminal_mismatch(for_seq, for_targ_seq, 
                                        three_prime_table, del_ct_thresh) or \
            is_fatal_terminal_mismatch(rev_seq, rev_targ_seq, 
                                        three_prime_table, del_ct_thresh):
                print("False Negative due to terminal MM failure!")
                print(assay_name)
                print(each_pos["Accession"])
                print(for_seq, for_targ_seq, rev_seq, rev_targ_seq)
                print()
                new_FN_list.append(each_pos)
        else:
            new_TP_list.append(each_pos)
        
    return ( assay_name, new_TP_list, new_FN_list)


def filter_three_prime(Full_Dict, three_prime_table, del_ct_thresh):

    print("\nRemoving 3' mismatches over threshold from True Positives...")
    assay_dict_list = Full_Dict["Assay"]
    for each_assay in assay_dict_list:
        FN_list = each_assay["False Negatives"]
        TP_list = each_assay["True Positives"]
        new_TP_list = []
        for each_pos in TP_list:

            ''' Get last two bases of Primers '''
            # Forward Primer
            for_prim = each_pos["Alignments"]["Forward Primer"]["5'"]
            for_seq = for_prim[-2:]

            # Reverse Primer
            rev_prim = each_pos["Alignments"]["Reverse Primer"]["5'"]
            rev_seq = rev_prim[-2:]

            ''' Get corresponding bases from Targets '''
            # Forward Primer
            for_prim_targ = each_pos["Alignments"]["Forward Primer"]["3'"]
            for_targ_seq = for_prim_targ[-2:]
            for_targ_seq = complement(for_targ_seq)

            # Reverse Primer
            rev_prim_targ = each_pos["Alignments"]["Reverse Primer"]["3'"]
            rev_targ_seq = rev_prim_targ[-2:]
            rev_targ_seq = complement(rev_targ_seq)

            ''' If terminal MM put in False Negative List; else put in new TP list '''
            if is_fatal_terminal_mismatch(for_seq, for_targ_seq, three_prime_table, del_ct_thresh) or \
                is_fatal_terminal_mismatch(rev_seq, rev_targ_seq, three_prime_table, del_ct_thresh):
                    print("False Negative due to terminal MM failure!")
                    print(each_assay["Name"])
                    print(each_pos["Accession"])
                    print(for_seq, for_targ_seq, rev_seq, rev_targ_seq)
                    print()
                    FN_list.append(each_pos)
            else:
                new_TP_list.append(each_pos)

        each_assay["True Positives"] = new_TP_list

    # Add timestamp
    Full_Dict.update({ "Timestamp" : str(datetime.now()) })

    return Full_Dict


def make_negatives_list(in_file, out_file):

    with open(in_file) as json_file:
        full_dict = json.load(json_file)

    with open(out_file, 'w') as write_file:
        for assay_dict in full_dict[assay]:
            this_assay = assay_dict[assay_name]
            FN_list = assay_dict[negatives_type]
            acc_list = []
            for negative in FN_list:
                seq_accesion = negative["Accession"]
                acc_list.append(seq_accesion)
            acc_string = "\t".join(acc_list)
            print("{}\t{}".format(this_assay, acc_string), file=write_file)


def make_negatives_list(in_file, out_file, fasta_dir):

    with open(in_file) as json_file:
        full_dict = json.load(json_file)

    assays_list = []
    for assay_dict in full_dict[assay]:
        this_assay = assay_dict[assay_name]
        FN_list = assay_dict[negatives_type]
        acc_list = []
        for negative in FN_list:
            seq_accesion = negative["Accession"]
            file = os.path.join(fasta_dir, "{}.fna".format(seq_accesion))

            with open(file, 'r') as read_file:
                header = read_file.readline().strip().strip('>')
            header_fields = header.split('|')
            isolate = header_fields[0]
            accession = header_fields[1]
            col_date = header_fields[2]
            if len(header_fields) >= 17:
                location = header_fields[9].split('/')[2]
                orig_lab = header_fields[14]
                sub_lab = header_fields[15]
                authors = header_fields[16]
            elif len(header_fields) >= 5:
                location = header_fields[4]
                orig_lab = ""
                sub_lab = ""
                authors = ""
            else:
                location = ""
                orig_lab = ""
                sub_lab = ""
                authors = ""
            metadata_list = [ col_date, location, orig_lab, sub_lab, authors ]
            acc_list.append({"Accession": seq_accesion, "Metadata": metadata_list})
        assay_fails = { "Assay": this_assay, "Accession_List": acc_list}
        assays_list.append(assay_fails)

    print("Writing {}...".format(out_file))
    with open(out_file, 'w') as write_file:
        header_line = "Assay\tAccession\tCollection_Date\tCountry\tOriginating_Lab\tSubmitting_Lab\tAuthors"
        print(header_line, file=write_file)
        for each_assay in assays_list:
            for each_acc in each_assay["Accession_List"]:
                meta = "\t".join(each_acc["Metadata"])
                line = "\t".join( [ each_assay[assay], each_acc["Accession"], meta ] )
                print(line, file=write_file)


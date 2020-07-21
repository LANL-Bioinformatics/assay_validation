#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 14:13:20 2020

@author: adanm
"""
import json
from datetime import datetime

def generate_db_stats(overlap_tallies, genbank_tallies, gisaid_tallies,
                      total_list_len, db_stats_file, db_totals_file):

    ''' ############ Prepare Dictionary of Stats for json file ############# '''
    filter_tally_dict = {}
    filter_tally_dict.update({"Overlap": {} })
    filter_tally_dict["Overlap"].update({"GenBank Initial": overlap_tallies[0]})
    filter_tally_dict["Overlap"].update({"GISAID Initial": overlap_tallies[1]})
    filter_tally_dict["Overlap"].update({"Exceptions": overlap_tallies[2]})
    filter_tally_dict["Overlap"].update({"Overlap": overlap_tallies[3]})
    filter_tally_dict["Overlap"].update({"GISAID Overlap Removed": overlap_tallies[4]})

    filter_tally_dict.update({"Filter GenBank": {} })
    filter_tally_dict["Filter GenBank"].update({"list before": genbank_tallies[0]})
    filter_tally_dict["Filter GenBank"].update({"removed": genbank_tallies[1]})
    filter_tally_dict["Filter GenBank"].update({"bat": genbank_tallies[2], "pangolin": genbank_tallies[3], "short": genbank_tallies[4]})
    filter_tally_dict["Filter GenBank"].update({"kept": genbank_tallies[5]})
    filter_tally_dict["Filter GenBank"].update({"list after": genbank_tallies[6]})

    filter_tally_dict.update({"Filter GISAID": {} })
    filter_tally_dict["Filter GISAID"].update({"list before": gisaid_tallies[0]})
    filter_tally_dict["Filter GISAID"].update({"removed": gisaid_tallies[1]})
    filter_tally_dict["Filter GISAID"].update({"bat": gisaid_tallies[2], "pangolin": gisaid_tallies[3], "short": gisaid_tallies[4]})
    filter_tally_dict["Filter GISAID"].update({"kept": gisaid_tallies[5]})
    filter_tally_dict["Filter GISAID"].update({"list after": gisaid_tallies[6]})

    filter_tally_dict.update({"Total list length": total_list_len})

    filter_tally_dict.update( {"Timestamp": str(datetime.now()) } )


    ''' ##################### Write db stats to File ####################### '''

    with open(db_stats_file, 'w') as write_file:
        print(json.dumps(filter_tally_dict, indent=4), file=write_file)



    ''' ############ Prepare Dictionary of db_totals.json file ############# '''
    db_total_dict = {
        "GISAID_tot": overlap_tallies[1],
        "GenBank_tot": overlap_tallies[0],
        "final_date": str(datetime.now()),
        "final_db_tot": total_list_len,
        "overlap": overlap_tallies[3]
    }

    ''' ##################### Write db_totals.json file ####################### '''

    with open(db_totals_file, 'w') as write_file:
        print(json.dumps(db_total_dict, indent=4), file=write_file)


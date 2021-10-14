#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 13:19:23 2021

@author: adanm
"""
from sys import platform
from datetime import datetime


def write_config_file(vars_ns, filenames_ns, txt):
    text_list = []
    
    ''' Write file header '''
    text_list.append(txt.startfile)
    text_list.append( datetime.now().strftime("Generated %A, %B %d, %Y: %H:%M:%S") )
    text_list.append( "\"\"\"\n" )

    '''' TNTBLAST location '''
    if platform == "darwin":
        text_list.append( "tntblast_location = \"/Users/adanm/.bin/\"")  # On Mac (spyder)
    elif platform == "linux":
        text_list.append( "tntblast_location = \"\"")     # On Rustang
        
    ''' Variables from Args '''
    text_list.append( "assay_type = '{}' # str, Assay format to simulate".format(vars_ns.assay_type) )
    text_list.append( "min_primer_Tm = {} # int, Minimum primer Tm".format(vars_ns.min_primer_Tm) )
    text_list.append( "min_probe_Tm = {} # int, Minimum probe Tm".format(vars_ns.min_probe_Tm) )
    text_list.append( "max_amplicon_length = {} # int, Maximum amplicon length".format(vars_ns.max_amplicon_length) )
    text_list.append( "minimum_target_length = {} ".format(vars_ns.minimum_target_length) +
                     "# int, Minimum length of target sequences. Doesn't filter GISAID" )
    text_list.append( "best_match = {} # Only".format(vars_ns.best_match) +
                     " save the best match, in Tm, between a query and target")
    text_list.append("\n\n")
    
    ''' Paths and Files '''
    text_list.append("# Paths")
    text_list.append("resource_directory =  \"{}\"".format(filenames_ns.resource_dir))
    text_list.append("results_directory =  \"{}\"".format(filenames_ns.results_dir))
    text_list.append("tnt_results_directory =  \"{}\"".format(filenames_ns.tnt_dir))
    text_list.append("fasta_directory =  \"{}\"\n".format(filenames_ns.fasta_dir))
    
    text_list.append("# Files")
    text_list.append("assay_file =  \"{}\"".format(filenames_ns.assay_file))
    text_list.append("sequence_file =  \"{}\"".format(filenames_ns.combined_fasta))

    ''' Write to file '''
    text_to_write = "\n".join(text_list)
    with open(filenames_ns.config_file, 'w') as write_file:
        print(text_to_write, file=write_file)
    return

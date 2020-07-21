#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:33:15 2020

@author: adanm
"""
import os
import re
import sys
import time
import argparse as ap


''' 'Global' Variables '''
# Data download switches
daily_gisaid_format = False
meta_require_format = True


''' Methods '''

def split_daily_gisaid_file(gisaid_file, download_dir):
    ''' ################################################ XXXX ### '''
    with open(gisaid_file, 'r') as file_handle:
        line = file_handle.readline()
        line2 = []
        fields = line[1:].split("/")
        subfields = fields[-1].split("|")
        filename = subfields[1] + ".fna"
        output_path = os.path.join(download_dir, filename)

        with open(output_path, "w") as write_handle:
            print(line.strip(), file=write_handle)


        for line in file_handle.readlines():
            if line[0] == ">":
                line2_str = "".join(line2)
                with open(output_path, "a") as write_handle:
                    print(line2_str.strip(), file=write_handle)
                line2 = []
                fields = line[1:].split("/")
                subfields = fields[-1].split("|")
                filename = subfields[1] + ".fna"
                output_path = os.path.join(download_dir, filename)
                with open(output_path, "w") as write_handle:
                    print(line.strip(), file=write_handle)
            else:
                line2.append(line.strip())


def get_header(isolate, metafile):

    header_dict = {}

    with open(metafile, 'r') as read_meta:
        column_headers = read_meta.readline().strip().split('\t')

        for metaline in read_meta.readlines():
            metaline_fields = metaline.strip().split('\t')
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



def split_gisaid_file_with_meta(gisaid_file, download_dir, metafile):
    ''' ################################################ XXXX ### '''
    with open(gisaid_file, 'r') as read_gisaid:
        line = read_gisaid.readline()
        isolate = line.strip().strip('>')
        line2 = []
        accession, header = get_header(isolate, metafile)

        if not accession == None:

            filename = accession + ".fna"
            output_path = os.path.join(download_dir, filename)


            with open(output_path, "w") as write_handle:
                print(header, file=write_handle)


        for line in read_gisaid.readlines():
            if line[0] == ">":
                line2_str = "".join(line2)
                with open(output_path, "a") as write_handle:
                    print(line2_str.strip(), file=write_handle)
                isolate = line.strip().strip('>')
                line2 = []
                accession, header = get_header(isolate, metafile)

                if not accession == None:

                    filename = accession + ".fna"
                    output_path = os.path.join(download_dir, filename)

                    with open(output_path, "w") as write_handle:
                        print(header, file=write_handle)
            else:
                line2.append(line.strip())


def get_arguments():
    parser = ap.ArgumentParser(prog=sys.argv[0].split('/')[-1],
                               description=sys.argv[0].split('/')[-1])
    parser.add_argument('-r', '--resource_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="download directory")
    parser.add_argument('-f', '--fasta_directory', metavar='[STR]', nargs=1, type=str,
                        required=True, help="download directory")
    return parser.parse_args()



def main():
    print("\nRunning:", sys.argv[0].split('/')[-1])
    tic = time.perf_counter()

    ''' Get command line argument '''
    args = get_arguments()

    resource_dir = args.resource_directory[0]
    fasta_dir = args.fasta_directory[0]


    # File Names
    metafile = os.path.join(resource_dir, "metadata.tsv")
    if daily_gisaid_format:
        gisaid_fasta_files = [ os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir)
                                if re.match(r'gisaid_hcov-19_2020_.*\.fasta$', f) ]
        gisaid_fasta_files.sort()
        gisaid_fasta_file = gisaid_fasta_files[-1]
    elif meta_require_format:
        gisaid_fasta_file = os.path.join(fasta_dir, "sequences.fasta")


    ''' ########################## GISAID Split ############################ '''


    ''' Split gisaid fasta into multiple fna files in target directory '''
    if os.path.exists(gisaid_fasta_file):
        print("\nSplitting GISAID sequence file...")
        if daily_gisaid_format:
            split_daily_gisaid_file(gisaid_fasta_file, fasta_dir)
        elif meta_require_format:
            split_gisaid_file_with_meta(gisaid_fasta_file, fasta_dir, metafile)
    else:
        print("\nNo Gisaid Fast File !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")


    toc = time.perf_counter()
    elapsed_time = toc - tic
    hours = int(elapsed_time / 3600)
    minutes = int(elapsed_time / 60 - hours * 60)
    seconds = elapsed_time - minutes * 60
    print("\nFinished")
    print("Total Elapsed Time: {} hours, {} minutes, {} seconds".format(hours, minutes, seconds))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import sys
import os
import re
import json

class AssayResult(object):
    """ This is the class for generating assay result JSON files.

    :param json: The file location of the input JSON tree
    :type json: str
    :param outdir: The output dir
    :type outdir: str
    """

    def __init__(self, ar_json, outdir):
        self.ar_json = ar_json
        self.outdir = outdir
        self.data = None
        self.stats = {
            "true_positives": {}
        }

    def load_json(self):        
        with open(self.ar_json) as f:
            self.data = json.load(f)

    def get_stats(self):        
        return self.stats

    def split_json(self, assay_list=None, genome_list=None):
            """ Split the whole assay result JSON file to "outdir/assay_id/genome_id.json".

            :param assay_list: Only create files for assay in this list
            :type assay_list: list
            :param genome_list: Only create files for genome in this list
            :type genome_list: list
            """

            # {
            #     "Assay": [
            #         {
            #             "False Negatives": [],
            #             "Name": "CDC-2019-nCoV_N1",
            #             "True Positives": [
            #                 {
            #                     "Accession": "MT253702",
            #                     "Alignments": {
            #                     ...

            for assay in self.data['Assay']:
                assay_id = assay['Name']

                # skip the assay if it's not in the assay_list
                if assay_list and (not assay_id in assay_list):
                    continue

                self.stats["true_positives"][assay_id] = {
                    'genome_num': 0,
                    'genome_source': {
                        'gisaid': 0,
                        'genbank': 0
                    },
                    'assay_sequence': {
                        "forward_primer": "",
                        "reverse_primer": "",
                        "probe": "",
                    }

                }

                if not os.path.exists(f'{assay_id}'):
                    os.makedirs(f'{self.outdir}/{assay_id}', exist_ok=True)

                genomes = assay["True Positives"]

                for genome in genomes:
                    genome_id = genome['Accession']

                    # stats
                    self.stats["true_positives"][assay_id]['genome_num'] += 1

                    if genome_id.startswith("EPI"):
                        self.stats["true_positives"][assay_id]['genome_source']['gisaid']+=1
                    else:
                        self.stats["true_positives"][assay_id]['genome_source']['genbank']+=1

                    # if not self.stats["true_positives"][assay_id]['assay_sequence']['forward_primer']:
                    #     self.stats["true_positives"][assay_id]['assay_sequence']['forward_primer'] = genome['Sequences']['Forward Primer']
                    #     self.stats["true_positives"][assay_id]['assay_sequence']['reverse_primer'] = genome['Sequences']['Reverse Primer']
                    #     self.stats["true_positives"][assay_id]['assay_sequence']['probe'] = genome['Sequences']['Probe']

                    # skip the genome if it's not going to display
                    if genome_id and (not genome_id in genome_list):
                        continue

                    with open(f'{self.outdir}/{assay_id}/{genome_id}.json', 'w') as f:
                        # print(f'writing {self.outdir}/{assay_id}/{genome_id}.json')
                        f.write(json.dumps(genome))
                    f.close()
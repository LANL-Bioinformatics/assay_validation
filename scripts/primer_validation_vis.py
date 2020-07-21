#!/usr/bin/env python3
import sys
import os
import argparse as ap
import json
import logging
from PhyHeatmap import PhyXML, Metadata
from AssayResult import AssayResult

__version__='0.0.7'

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M',
)

def parse_params():
    p = ap.ArgumentParser(prog='primer_validation_vis.py',
                          description="""Generate XML for visualization of primer validation""")

    p.add_argument('-n', '--newick',
                   metavar='[FILE]', type=str, required=True,
                   help="Input newick tree")

    p.add_argument('-o', '--outfile',
                   metavar='[XML]', type=str, required=True,
                   help="Output XML file")

    p.add_argument('-p', '--primer',
                   metavar='[FILE]', type=str, required=True,
                   help="Input validation file for primers")

    p.add_argument('-a', '--assayseq',
                   metavar='[FILE]', type=str, required=True,
                   help="Assay sequence file")

    p.add_argument('-r', '--resultjson',
                   metavar='[FILE]', type=str, required=False,
                   help="Assay results in JSON format")

    p.add_argument('-m', '--metadata',
                   metavar='[JSON]', nargs='+', type=str, required=True,
                   help="Input metadata JSON for all strains")

    args_parsed = p.parse_args()
    return args_parsed


def main():
    argvs = parse_params()

    logging.info("Loading metadata...")
    metadata = Metadata(argvs.metadata)

    logging.info("Initing tree...")
    phyxml = PhyXML(argvs.newick,
                    argvs.outfile,
                    argvs.assayseq,
                    metadata,
                    argvs.primer)

    logging.info("Removing leaves not in the match table...")
    phyxml.purge_tree()

    logging.info("Collapsing leaves by branch leangth and patterns...")
    phyxml.collapsing_nodes(collapse_node_by_branch=True,
                            collapse_node_by_pattern=True,
                            collapse_branch_len=0)

    logging.info("Adding metadata to nodes...")
    phyxml.prepare_ext_phyloxml()

    logging.info("Processing validation results...")
    g_dom_list = []

    # prepare graphics dom
    # title = os.path.basename(argvs.primer)
    g_dom = phyxml.add_graph_validation_heatmap(argvs.primer, "Validation result")
    g_dom_list.append(g_dom)

    logging.info("Adding graphics...")
    # g_dom = phyxml.add_graph_bars("Stats")
    # phyxml.append_dom_to_tag('graphs', g_dom)

    for dom in g_dom_list:
        phyxml.append_dom_to_tag('graphs', dom)

    # adding colors to taxonomies
    ts_dom = phyxml.add_taxanomy_color()
    phyxml.append_dom_to_tag('phyloxml', ts_dom)

    logging.info("Writing XML...")
    phyxml.writeXML()

    logging.info("Writing metadata JSON...")
    phyxml.writeMetaJson()

    logging.debug("Continental colors...")
    logging.debug(phyxml.cont_colors)

    # Parsing assay result JSON
    outdir = os.path.dirname(argvs.outfile)
    ar = AssayResult(argvs.resultjson, f'{outdir}/assay_result_json')

    logging.info("Loading assay results JSON...")
    ar.load_json()

    logging.info("Spliting JSON...")
    assay_list = phyxml.get_assay_list()
    genome_list = phyxml.get_genome_list()
    ar.split_json(assay_list, genome_list)

    logging.info("Writing stats JSON...")
    from datetime import datetime
    stats = {
        "time": str(datetime.now()),
        "tree": phyxml.get_stats(),
        "validation": ar.get_stats()
    }
    with open(f'{argvs.outfile}.stats.json', 'w') as f:
        f.write(json.dumps(stats))

    logging.info("Completed.")


if __name__ == "__main__":
    logging.info(f'Assay Validation Visualization v{__version__}')
    main()
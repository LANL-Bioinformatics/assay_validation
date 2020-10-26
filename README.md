# COVID-19 In silico Evaluation of Diagnostic Assays

This is a web-app and workflow for computationally screening PCR-based assays against the increasing number of SARS-CoV-2 sequences being deposited for public release in GISAID and GenBank. Assays are assessed for binding based on free energy and melting temperature to determine whether binding occurs between the assay oligonucleotides (primers and probe), and target sequence.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

---

**NOTE:** 
The program is designed to retrieve genomes from GISAID and/or NCBI for the users' convenience. We have removed any GISAID data from our web-app per GISAID's request since September 2020.

---

## LANL software requirements

The following software are required and the executable binaries should be put under your system PATH:

1. ThermonucleotideBLAST
    
    ThermonucleotideBLAST is a software program for searching a target database of nucleic acid sequences using an assay-specific query. The detail instructions for installation can be found on software's github https://github.com/jgans/thermonucleotideBLAST.

2. EpiCoV_downloader

    This is a GISAID sequence downloader for retrieving SARS-COV-2 genomes and the metadata. The detail instructions for installation can be found on software's github https://github.com/poeli/EpiCoV_downloader.

3. PhyD3-am

    This web-app includes a modified version of PhyD3 (GPL3) and visualizations specifically developed for displaying stats, metadata, phylogenetic tree and assay evaluation results. Here is a simple installation and usage:

    ```
    $ cd phyd3-am/
    $ npm install
    $ node phyd3.js
    ```

    The detail information can be found on original software's github https://github.com/vibbits/phyd3. 

## Python library requirements

The following Python libraries are required:

* biopython >= 1.0
* bokeh  >= 1.0
* ete3 >= 3.1
* pandas >= 0.24.0


## Additional File Requirements

These files must be present in the resource directory (Indicated by the `-r` argument):
* "assays.txt"

    A list of assays with corresponding oligo sequences.  Each to a line.
Order of Oligo sequences: Forward primer Reverse primer, Probe.  Ex:

```
CDC-2019-nCoV_N1 GACCCCAAAATCAGCGAAAT TCTGGTTACTGCCAGTTGAATCTG ACCCCGCATTACGTTTGGTGGACC
```
* "reduced_assays.txt"

    A reduced list of assays to be used in the final results.  Same format
as "assays.txt"
* "del_ct_table.txt"

    A table of delta Ct values from { Li B, Kadura I, Fu DJ, Watson DE.
Genotyping with TaqMAMA. Genomics. 2004 Feb 1;83(2):311-20. }.
Tab separated.  Headers for Rows and columns.  First entry: "Row".  
Example first two lines:

```
Row     CC      GC      AC      TC      CG      GG      AG      TG      CA      GA      AA      TA      CT      GT      AT      TT
CC      0.0     0.3     0.5     -0.3    6.3     17.5    19.2    11.9    5.1     15.7    11.4    12.4    0.4     11.0    10.4    3.7
```


## Workflow

The workflow mainly includes 3 major steps. All scripts mentioned below can be found in the `scripts/` directory.

1. Download genomes from the two repositories and be cross-validated to remove any duplicate entries.

    The `am_download.py` script will download and prepare the necessary data for running the `assay_monitor.py` script.  Other scripts are called which log into the GISAID database to download SARS-CoV-2 genomes as individual fasta (.fna) files. NCBI is also accessed using NCBI's e-Utilities to download fastas from GenBank. Metadata for both GISAID and GenBank are downloaded from these sources.

2. Assess assays for binding based on free energy and melting temperature to determine whether binding occurs between the assay oligonucleotides (primers and probe), and target sequence.

    The `assay_monitor.py` runs the assay monitor functionality for this package. Using the fasta files downloaded by `am_download.py` and certain resource files (detailed below), assays (provided in a resource file) will be evaluated against each genome (previously downloaded fastas) using ThermonucleotideBLAST. Results will be output as full results in `Assay_Results.json`, as a summary in `summary_table.json`, and as cross referenced results in `match_table.csv`. Summary stats are also output to `db_stats.json` and `db_totals.json`. These files are used downstream in this package for producing visualizations.
   
3. Integrate the input phylogenetic tree with the assay evaluation results, then generate essential files for visualization.

    This tree visualization is rendered using a custom PhyD3 phylogenetic tree viewer. The `primer_validation_vis.py` will take an phylogenetic tree (we recommand PhaME) and the output generated from `assay_monitor.py` script to produce essential data for the web-app to `<PhyD3_path>/dist/data/`. The data include heatmap-associated tree in extended phyloXML format along with results and stats of assay evaluation in individual JSON files for user's information. The heatmap displays the predicted mismatches and assay outputs for each genome of SARS-CoV-2. 

## Usage

We provide a script `am_start.sh` to glue all the scripts and to copy essential files to the web-app directory. Please use `am_start.sh -h` for more details.

`usage: am_start.sh [-h] -u [STR] -p [STR] -e [STR] -m [STR] -M [STR] -s [STR] -r [STR] -R [STR] -f [STR] -t [STR]`

## Citation

Li, P.E., Myers y Gutiérrez, A., Davenport, K., Flynn, M., Hu, B., Lo, C.C., Jackson, E.P., Shakya, M., Xu, Y., Gans, J. and Chain, P.S., 2020. A Public Website for the Automated Assessment and Validation of SARS-CoV-2 Diagnostic PCR Assays. arXiv preprint arXiv:2006.04566.

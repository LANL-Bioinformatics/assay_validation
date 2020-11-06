#!/bin/bash

# This is a script to run assay monitor workflow
# Po-E Li, B-10, LANL

##### Functions
init_dir()
{
    local DIRNAME=$1
    [ ! -d $DIRNAME ] && mkdir -p $DIRNAME;
    [ ! -d $DIRNAME ] && { echo "ERROR: Directory $DIRNAME DOES NOT exists."; exit 1; }
}

check_empty_argument()
{
    local ARG=$1
    local VAL=$2
    [ -z $VAL ] && { echo "ERROR: --$ARG not specified"; usage; exit 1; }
}

usage()
{
    echo "usage: am_start.sh [-h] -e [STR] -m [STR] -M [STR] -s [STR] -r [STR] -R [STR] -f [STR] -t [STR] [ -D ]

The following arguments are required:

    '-e', '--email' <email_address>'
        This is a requirement from ncbi for inclusion in a query

    '-r', '--resource_directory' <resource_directory>
        This is the directory to which certain resource files will be
        created.
        Use the format:  /XXX/XXX/XXX/

    '-R', '--results_directory' <results_directory>
        This is the directory to which certain results files will be
        downloaded or created.
        Use the format:  /XXX/XXX/XXX/

    '-f', '--fasta_directory' <download_directory>
        This is the directory to which fna sequence files will be
        downloaded.
        Use the format:  /XXX/XXX/XXX/

    '-t', '--tnt_results_directory' <tnt_results_directory>
        This is the directory to which results files will be created when
        TNTBLAST is run.
        Use the format:  /XXX/XXX/XXX/

    '-d', '--phyd3d_dist_directory' <dist_directory>
        This is the path of dist/ directory under the PhyD3 installation

    '-m', '--mindate' <early_date>
        Put the earliest publication date or lower range of date in format:
        'YYYY/MM/DD'

    '-M', '--maxdate' <late_date>
        Put the latest publication date or upper range of date in format:
        'YYYY/MM/DD'

    '-s', '--seqlengthrange' <seq_length_range>
        Input the range of sequence lengths in the format:  29000:30500.
        For no limit in either direction, put: 0:99999999

    '-n', '--newick_treefile' <newick_tree_file>
        Input a tree file in newick format


The following arguments are optional:
    '-D', '--dont_download'
        Use this flag to skip downloading sequences.  Sequences that 
        have already been downloaded to the fasta_directory will be
        processed by the assay validation script.


"
}

##### Main
email=
mindate=
maxdate=
seqlengthrange=
newick_treefile=
resource_directory=
results_directory=
fasta_directory=
tnt_results_directory=
phyd3d_dist_directory='phyd3-am/dist'
skip_download=false

while [ "$1" != "" ]; do
    case $1 in
        -e | --email )              shift
                                    email=$1
                                    ;;
        -m | --mindate )            shift
                                    mindate=$1
                                    ;;
        -M | --maxdate )            shift
                                    maxdate=$1
                                    ;;
        -s | --seqlengthrange )     shift
                                    seqlengthrange=$1
                                    ;;
        -n | --newick_treefile )    shift
                                    newick_treefile=$1
                                    ;;
        -r | --resource_directory ) shift
                                    resource_directory=$1
                                    ;;
        -R | --results_directory )  shift
                                    results_directory=$1
                                    ;;
        -f | --fasta_directory )    shift
                                    fasta_directory=$1
                                    ;;
        -t | --tnt_results_directory )  shift
                                        tnt_results_directory=$1
                                        ;;
        -d | --phyd3d_dist_directory )  shift
                                        phyd3d_dist_directory=$1
                                        ;;
        -D | --dont_download )      skip_download=true
                                    ;;
        -h | --help )               usage
                                    exit
                                    ;;
        * )                         usage
                                    exit 1
    esac
    shift
done

# verify arguments
check_empty_argument email $email
check_empty_argument mindate $mindate
check_empty_argument maxdate $maxdate
check_empty_argument seqlengthrange $seqlengthrange
check_empty_argument newick_treefile $newick_treefile
init_dir $resource_directory
init_dir $results_directory
init_dir $fasta_directory
init_dir $tnt_results_directory
init_dir $phyd3d_dist_directory/data

set -e;
export PATH=scripts/:$PATH

if ! $skip_download; then
    echo -e "\nDownloading genomes...\n"
    am_download.py \
        -e $email \
        -m $mindate \
        -M $maxdate \
        -s $seqlengthrange \
        -r $resource_directory \
        -R $results_directory \
        -f $fasta_directory
else
    echo -e "\nSkipping download\n"
fi


echo -e "\nEvaluating assays...\n"
assay_monitor.py \
    -r $resource_directory \
    -R $results_directory \
    -f $fasta_directory \
    -t $tnt_results_directory

# generate genbank metadata
echo "isolate|accession|col_date|create_date|Country: region|seq_length|start_pos|end_pos|region|country|division|location, region_exposure|country_exposure|division_exposure|segment|host|age|sex|originating_lab|submitting_lab|authors|title|comment" | sed "s/|/\t/g" > $resource_directory/metadata.tsv
cat $fasta_directory/All_Seqs.fasta | grep '>' | sed "s/>//" | sed "s/|/\t/g" >> $resource_directory/metadata.tsv

echo -e "\nGenerating data for visualization...\n"
primer_validation_vis.py \
    -n $newick_treefile \
    -m $resource_directory/metadata.tsv \
    -p $results_directory/match_table.csv \
    -a $resource_directory/assays.txt \
    -r $results_directory/Assay_Results.json \
    -o $phyd3d_dist_directory/data/SARS-CoV-2.xml

cp $resource_directory/db_totals.json $phyd3d_dist_directory/data/
cp $results_directory/summary_table.json $phyd3d_dist_directory/data/

set +xe;

echo -e "\nDone."
echo "Open http://localhost:8080/ at your browser! Enjoy!"

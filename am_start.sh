#!/bin/bash

# This is a script to run assay monitor workflow
# Po-E Li, B-10, LANL

##### Functions
init_dir()
{
    local DIRNAME=$1
    [ -d $DIRNAME ] && mkdir -p $DIRNAME;
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
    echo "usage: am_start.sh [-h] -u [STR] -p [STR] -e [STR] -m [STR] -M [STR] -s [STR] -r [STR] -R [STR] -f [STR] -t [STR]

The following arguments are required: 
    '-u', '--gisaid_username' <username>
        This is your username to log into gisaid.org

    '-p', '--gisaid_passwd' <password>
        This is your password for your GISAID login

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

"
}

##### Main
gisaid_username=
gisaid_passwd=
email=
mindate=
maxdate=
seqlengthrange=
newick_treefile=
resource_directory=
results_directory=
fasta_directory=
tnt_results_directory=
phyd3d_dist_directory=

while [ "$1" != "" ]; do
    case $1 in
        -u | --gisaid_username )    shift
                                    gisaid_username=$1
                                    ;;
        -p | --gisaid_passwd )      shift
                                    gisaid_passwd=$1
                                    ;;
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
        -h | --help )               usage
                                    exit
                                    ;;
        * )                         usage
                                    exit 1
    esac
    shift
done

# verify arguments
check_empty_argument gisaid_username $gisaid_username
check_empty_argument gisaid_passwd $gisaid_passwd
check_empty_argument email $email
check_empty_argument mindate $mindate
check_empty_argument maxdate $maxdate
check_empty_argument seqlengthrange $seqlengthrange
check_empty_argument newick_treefile $newick_treefile
init_dir $resource_directory
init_dir $results_directory
init_dir $fasta_directory
init_dir $tnt_results_directory
init_dir $phyd3d_dist_directory

echo "Downloading genomes..."
am_download.py \
    -u $gisaid_username \
    -p $gisaid_passwd \
    -e $email \
    -m $mindate \
    -M $maxdate \
    -s $seqlengthrange \
    -r $resource_directory \
    -R $results_directory \
    -f $fasta_directory \
    -t $tnt_results_directory

echo "Evaluating assays..."
assay_monitor.py \
    -r $resource_directory \
    -R $results_directory \
    -f $fasta_directory \
    -t $tnt_results_directory

echo "Generating data for visualization..."
primer_validation_vis.py \
    -n $newick_treefile \
    -m $resource_directory/metadata.tsv \
    -p match_table.csv \
    -a $resource_directory/assay.txt \
    -r Assay_Results.json \
    -o dist/data/SARS-CoV-2.xml

cp $resource_directory/db_totals.json dist/data/
cp $resource_directory/summary_table.json dist/data/

echo "Done."
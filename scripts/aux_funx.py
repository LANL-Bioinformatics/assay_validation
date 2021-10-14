import re

def get_vname_generator_from(sequence_file, iso_idx):
    if iso_idx == 0: # case gisaid
        start_idx = 1 # trim ">" from common name of isolate
    else:            # case genbank
        start_idx = 0 # don't trim accession number

    with open(sequence_file, 'r') as read_file:
        gis_lines = read_file.readlines()

    return (
        line.strip().split("|")[iso_idx][start_idx:]
        for line in gis_lines
        if line[0] == ">"
    )


def virus_acc_map_from(metafile):

    # String methods needed by this method: #
    def virus_name(metaline, idx):
        metaline_fields = metaline.strip().split('\t')
        return metaline_fields[idx]

    def virus_name_spaceless(metaline, idx):
        metaline_fields = metaline.strip().split('\t')
        return re.sub(r"\s+", "", metaline_fields[idx])

    def gisaid_accession(metaline, idx):
        metaline_fields = metaline.strip().split('\t')
        return metaline_fields[idx]

    # Begin this method: #
    with open(metafile, 'r') as read_meta:
        # Read column header and get indices for V_name and accession
        column_headers = read_meta.readline().strip().split('\t')
        v_name_idx = column_headers.index("Virus name")
        acc_idx = column_headers.index("Accession ID")

        # Load remainder of file
        meta_lines = read_meta.readlines()

    # Generate virus_name: accession mapping
    virus_acc_map = {
        virus_name(metaline, v_name_idx): gisaid_accession(metaline, acc_idx)
        for metaline in meta_lines
    }
    # Generate spaceless virus_name: accession mapping
    spaceless_virus_acc_map = {
        virus_name_spaceless(metaline, v_name_idx
            ): gisaid_accession(metaline, acc_idx)
        for metaline in meta_lines
    }
    # Join dicts
    virus_acc_map.update(spaceless_virus_acc_map)
   
    return virus_acc_map


def create_accession_list(vname_generator, vname_acc_map):
    return [
        vname_acc_map[vname]
        for vname in vname_generator
    ]


def write_acc_list(accession_list, accession_file):
    with open(accession_file, 'w') as write_file:
        for accession in accession_list:
            print(accession, file=write_file)


def generate_accession_file(file_ns, iso_idx):
    # Get list of virus names from gisaid file
    vname_generator = get_vname_generator_from(file_ns.sequence_file, iso_idx)
    if iso_idx == 0: # Do this only for gisaid metadata
        # Map virus names to accessions using metadata
        vname_acc_map = virus_acc_map_from(file_ns.metafile)
        # Create accession list and write to file
        accession_list = create_accession_list(vname_generator, vname_acc_map)
    else:
        accession_list = list(vname_generator)
    write_acc_list(accession_list, file_ns.accession_file)
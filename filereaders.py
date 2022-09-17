


def readInputFile(input_file_path):

    with open(input_file_path) as f:
        lines = f.readlines()

    clean_lines = [line.strip() for line in lines]

    return clean_lines

def readPanelFile(input_file_path):

    clean_lines = readInputFile(input_file_path)
    panel={}

    for line in clean_lines:

        spat = line.split(":")
        primer_set=spat[0]
        loci_in_primer_set=spat[1].split(",")
        loci_info = {}
        for each_loci in  loci_in_primer_set:

            spat = each_loci.split("(")
            msi_loci_name=spat[0].strip()
            msi_loci_start = spat[1].strip()
            #msi_loci_end = spat[2].strip()

            loci_info[msi_loci_name] = [msi_loci_start,0]# msi_loci_end]

        panel[primer_set] = loci_info

    return panel

def figure_out_loci_from_run_name(panel, run_name):

    for primer_set in panel.keys():
        if primer_set in run_name:
            return primer_set

    return "FAIL"
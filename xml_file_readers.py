import xml.etree.ElementTree as ET


def readPanelXml(input_file_path):

    mytree = ET.parse(input_file_path)
    myroot = mytree.getroot()

    primer_sets_data_by_primer_set_name={}

    for primer_set in myroot:
        primer_set_name=primer_set.attrib["name"]
        primer_set_data_by_loci={}
        for loci in primer_set:
            loci_name = loci.attrib["name"]
            loci_data={}

            for data in loci:

                incoming_txt=data.text.strip()
                incoming_tag = data.tag.strip()

                if (incoming_tag == "length"):
                    interval_splat = incoming_txt.split('-')
                    msi_loci_start = int(interval_splat[0].strip())
                    msi_loci_end = int(interval_splat[1].strip())
                    loci_data["length"] = [msi_loci_start, msi_loci_end]
                else:
                    loci_data[incoming_tag]=incoming_txt

            primer_set_data_by_loci[loci_name]=loci_data

        primer_sets_data_by_primer_set_name[primer_set_name]=primer_set_data_by_loci

    print("Panel Data Loaded:")
    print(primer_sets_data_by_primer_set_name)

    return primer_sets_data_by_primer_set_name

def figure_out_loci_from_run_name(panel, run_name):

    for primer_set in panel.keys():
        if primer_set in run_name:
            return panel[primer_set]

    return False


def read_truth_data(input_file_path):

    mytree = ET.parse(input_file_path)
    myroot = mytree.getroot()

    truthdata_by_sample_name={}

    for sample in myroot:

        sample_name=sample.attrib["name"]
        sample_data_by_loci = {}
        for loci in sample:
            loci_name = loci.attrib["name"]
            loci_data= []

            for data in loci:

                incoming_txt=data.text.strip()
                incoming_tag = data.tag.strip()

                if (incoming_tag == "alleles"):
                    alleles_splat = incoming_txt.split(',')
                    loci_data = [int(a) for a in alleles_splat]

            sample_data_by_loci[loci_name]=loci_data

        truthdata_by_sample_name[sample_name]=sample_data_by_loci

    print("Truth Data Loaded:")
    print(truthdata_by_sample_name)

    return truthdata_by_sample_name

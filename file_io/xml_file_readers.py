import os
import xml.etree.ElementTree as ET

class TruthDataForSample:

    def __init__(self, species_name, sample_data_by_loci):
        self.truth_by_loci = sample_data_by_loci
        self.species_name = species_name

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

def figure_out_loci_from_file_name(panel, file_name):

    base = os.path.basename(file_name).split(".")[0]
    primer_set_equivalence={}
    primer_set_equivalence["PS1"]="I3_B20_A1"
    primer_set_equivalence["PS2"]="B11_I14_C8"
    primer_set_equivalence["PS3"]="B9_B18_E9"
    primer_set_equivalence["PS4"]="BF3_BF19_B6"
    primer_set_equivalence["PS5"]="BF15_B266_A3"
    primer_set_equivalence["PSM"]="DoesntExistInMWData"
    primer_set_equivalence["PSE"]="DoesntExistInMWData"

    for primer_set in panel:
        if primer_set in base :
            return panel[primer_set]
        if primer_set_equivalence[primer_set] in base :
            return panel[primer_set]


    return False


def read_truth_data(input_file_path):

    mytree = ET.parse(input_file_path)
    myroot = mytree.getroot()

    truthdata_by_sample_name={}

    for sample in myroot:

        sample_name = sample.attrib["name"]

        if 'species' in sample.attrib:
            species_name = sample.attrib['species']
        else:
            species_name = False

        sample_data_by_loci = {}
        for loci in sample:
            loci_name = loci.attrib["name"]
            loci_data= []
            #print(loci_name )
            for data in loci:

                #print(data.text)
                incoming_txt=data.text.strip()
                incoming_tag = data.tag.strip()

                if (incoming_tag == "alleles"):
                    alleles_splat = incoming_txt.split(',')
                    loci_data = [int(a) for a in alleles_splat]

            sample_data_by_loci[loci_name]=loci_data

        truthdata_by_sample_name[sample_name]=TruthDataForSample(species_name,sample_data_by_loci)

    print("Truth Data Loaded:")
    print(truthdata_by_sample_name)

    return truthdata_by_sample_name

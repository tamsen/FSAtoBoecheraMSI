import os
import accuracy
from datetime import datetime

import log


def write_results(outputDir, data):
    resultsFile = os.path.join(outputDir, "ResultsByFile.txt")
    now = datetime.now()
    day = now.strftime("%d/%m/%Y")
    time = now.strftime("%H:%M:%S")
    data_string = ",".join(data)
    time_stamp_string = ",".join([day, time])

    with open(resultsFile, 'a') as f:
        f.write(time_stamp_string + "," + data_string + "\n")


def consolidate_by_file_results_to_by_sample_results(results_by_file, panel_info, truth_info):
    FSA_results_by_sample_by_loci = {}
    primer_sets = panel_info.keys()

    for file in results_by_file.keys():

        sample_name = get_sample_name_from_file_name(file, primer_sets)

        if not (sample_name in FSA_results_by_sample_by_loci):
            FSA_results_by_sample_by_loci[sample_name] = {}

        truth_for_this_sample = accuracy.find_truth_for_this_sample(sample_name, truth_info)

        if not results_by_file[file]:
            continue

        for loci in results_by_file[file].MSI_loci_results_by_loci.keys():

            msi_results_for_loci = results_by_file[file].MSI_loci_results_by_loci[loci]

            if truth_for_this_sample:

                truth_by_loci = truth_for_this_sample.truth_by_loci
                true_species = truth_for_this_sample.species_name

                if truth_by_loci and loci in truth_by_loci:
                    truth_for_loci = truth_by_loci [loci]
                    if len(truth_for_loci) > 0:
                        msi_results_for_loci.set_truth_and_accuracy(truth_for_loci,true_species)

            # Note: special handling for PSE here, where we had to run a second FSA file
            # to rescue the E9 loci that failed in the original PS3
            loci_already_processed_for_this_sample=FSA_results_by_sample_by_loci[sample_name].keys()
            if loci == "E9":
                if("PS3" in file) and ("E9" in loci_already_processed_for_this_sample):
                    log.write_to_log("Processing file " + file )
                    log.write_to_log("Loci already called: " + str(loci_already_processed_for_this_sample))
                    prior_E9_data = FSA_results_by_sample_by_loci[sample_name]["E9"].final_alleles_called
                    log.write_to_log("Prior E9 calls: " + str(prior_E9_data))
                    log.write_to_log("Will not override prior E9 data. " )
                else:
                    FSA_results_by_sample_by_loci[sample_name][loci] = msi_results_for_loci
            else:
                FSA_results_by_sample_by_loci[sample_name][loci] = msi_results_for_loci

    return FSA_results_by_sample_by_loci


def get_sample_name_from_file_name(file, primer_sets):
    base = os.path.basename(file).split(".")[0].split("_")[0]
    for primer_set in primer_sets:
        base = base.replace(primer_set, "")
    base = base.split("-")[0]
    sample_name = base
    return sample_name


def write_summary_file(outputDir, bySampleResults, panel_info, final_calls):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])

    if final_calls:
        summaryFile = os.path.join(outputDir, "FinalCallsBySample" + time_stamp_string + ".tsv")
    else:
        summaryFile = os.path.join(outputDir, "RawCallsBySample" + time_stamp_string + ".tsv")

    primer_sets = panel_info.keys()
    expected_space_for_calls = 5

    header1_data = ["PrimerSet->"]
    header2_data = ["Loci->"]
    for primer_set in primer_sets:
        for loci in panel_info[primer_set].keys():
            for i in range(0, expected_space_for_calls):
                header1_data.append(primer_set)
                header2_data.append(loci)

    header1 = "\t".join([str(p) for p in header1_data])
    header2 = "\t".join([str(p) for p in header2_data])

    with open(summaryFile, 'w') as f:

        f.write(header1 + "\n")
        f.write(header2 + "\n")

        # for file in results_by_file:
        for sample in bySampleResults.keys():

            data_list = [sample]
            results_for_file = bySampleResults[sample]
            if (len(results_for_file.keys()) < 1):
                data_list.append("Analysis fail. No alleles detected")
            else:

                for primer_set in primer_sets:
                    for loci in panel_info[primer_set].keys():

                        if loci in results_for_file:

                            if final_calls:
                                allele_calls = results_for_file[loci].final_alleles_called
                            else:
                                allele_calls = results_for_file[loci].raw_alleles_called
                        else:
                            allele_calls = ["" for x in range(0, expected_space_for_calls)]

                        for i in range(0, expected_space_for_calls):

                            if (i < len(allele_calls)):
                                data_list.append(str(allele_calls[i]))

                            else:
                                data_list.append("-")

            data_line = "\t".join(data_list)
            f.write(data_line + "\n")

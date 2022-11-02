import os
from datetime import datetime


def assess_accuracy(outputDir, results_by_file, panel_info, truth_info):

        now = datetime.now()
        day = now.strftime("%d_%m_%Y")
        time = now.strftime("%H_%M_%S")
        time_stamp_string = "_".join([day, time])
        summaryFile = os.path.join(outputDir, "Accuracy_" + time_stamp_string + ".tsv")

        primer_sets = panel_info.keys()
        #expected_space_for_calls = 5

        header1_data = ["PrimerSet->"]
        header2_data = ["Loci->"]
        for primer_set in primer_sets:
            for loci in panel_info[primer_set].keys():
                #for i in range(0, expected_space_for_calls):
                    header1_data.append(primer_set)
                    header2_data.append("observed " + loci)
                    header2_data.append("expected " + loci)

        header1 = "\t".join([str(p) for p in header1_data])
        header2 = "\t".join([str(p) for p in header2_data])

        #bySampleResults = consolite_by_file_results_to_by_sample_results(results_by_file, panel_info)

        with open(summaryFile, 'w') as f:

            f.write(header1 + "\n")
            f.write(header2 + "\n")

            # for file in results_by_file:
            for file_path in results_by_file.keys():

                have_truth_for_this_sample = False

                # Get filename from FSA file path
                file_name = os.path.basename(file_path)
                # Figure out if we have truth data available for that sample

                for sample in truth_info.keys():
                    if sample in file_name:
                        truth_for_this_sample = find_truth_for_this_sample(sample, truth_info)


                data_list = [file_path]
                results_for_file = results_by_file[file_path]
                if len(results_for_file.keys()) < 1:
                    data_list.append("Analysis fail. No alleles detected")
                elif have_truth_for_this_sample  == False:
                    data_list.append("Can't determine accuracy. No truth data for this sample.")
                else:

                    for primer_set in primer_sets:
                        for loci in panel_info[primer_set].keys():

                            expected_alleles = truth_for_this_sample[loci]

                            if loci in results_for_file:
                                allele_calls = results_for_file[loci]
                            else:
                                #allele_calls = ["" for x in range(0, expected_space_for_calls)]
                                allele_calls = ["-"]

                            data_list.append( str(allele_calls))
                            data_list.append( str(expected_alleles))
                            #data_list.append(str(allele_calls) + ',' + str(expected_alleles))

                data_line = "\t".join(data_list) + "\n"
                f.writelines([data_line])


def find_truth_for_this_sample(sample, truth_info):

    truth_for_this_sample= False

    if sample in truth_info.keys():
        truth_for_this_sample = truth_info[sample]
        return truth_for_this_sample

    second_try= sample.split("c")[0]
    if second_try in truth_info.keys():
        truth_for_this_sample = truth_info[second_try]
        return truth_for_this_sample

    third_try= sample.split("C")[0]
    if third_try in truth_info.keys():
        truth_for_this_sample = truth_info[third_try]

    return truth_for_this_sample

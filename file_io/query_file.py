import os
import accuracy
from datetime import datetime

import log


def write_query_file(outputDir, bySampleResults, panel_info, final_calls):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])

    if final_calls:
        summaryFile = os.path.join(outputDir, "BatchQuery_FinalCallsBySample" + time_stamp_string + ".tsv")
    else:
        summaryFile = os.path.join(outputDir, "BatchQuery_RawCallsBySample" + time_stamp_string + ".tsv")

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

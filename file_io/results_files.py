import os
import accuracy
from datetime import datetime

def write_results(outputDir, data):

        resultsFile = os.path.join(outputDir,"ResultsByFile.txt")
        now = datetime.now()
        day = now.strftime("%d/%m/%Y")
        time = now.strftime("%H:%M:%S")
        data_string = ",".join(data)
        time_stamp_string = ",".join([day,time])

        with open(resultsFile, 'a') as f:
           f.write(time_stamp_string + "," + data_string + "\n")

def consolidate_by_file_results_to_by_sample_results(results_by_file, panel_info, truth_info):

        FSA_results_by_sample_by_loci={}
        primer_sets = panel_info.keys()

        for file in results_by_file.keys():

                base = os.path.basename(file).split(".")[0].split("_")[0]
                for primer_set in primer_sets:
                        base = base.replace(primer_set,"")

                base = base.split("--")[0]
                sample_name=base

                if not(sample_name in FSA_results_by_sample_by_loci):
                        FSA_results_by_sample_by_loci[sample_name] = {}

                truth_for_this_sample = accuracy.find_truth_for_this_sample(sample_name, truth_info)

                for loci in results_by_file[file].MSI_loci_results_by_loci.keys():

                        truth_for_loci = truth_for_this_sample[loci]
                        msi_results_for_loci=results_by_file[file].MSI_loci_results_by_loci[loci]
                        msi_results_for_loci.truth_data = truth_for_loci

                        FSA_results_by_sample_by_loci[sample_name][loci] = msi_results_for_loci


        return FSA_results_by_sample_by_loci

def write_summary_file(outputDir, bySampleResults , panel_info):

        now = datetime.now()
        day = now.strftime("%d_%m_%Y")
        time = now.strftime("%H_%M_%S")
        time_stamp_string = "_".join([day,time])
        summaryFile = os.path.join(outputDir,"ResultsBySample" + time_stamp_string +".tsv")

        primer_sets=panel_info.keys()
        expected_space_for_calls = 5

        header1_data=["PrimerSet->"]
        header2_data=["Loci->"]
        for primer_set in primer_sets:
                for loci in panel_info[primer_set].keys():
                        for i in range(0,expected_space_for_calls ):
                                header1_data.append(primer_set)
                                header2_data.append(loci)

        header1 = "\t".join([str(p) for p in header1_data])
        header2 = "\t".join([str(p) for p in header2_data])

        with open(summaryFile , 'w') as f:

                f.write(header1 + "\n")
                f.write(header2 + "\n")

                #for file in results_by_file:
                for sample in bySampleResults.keys():

                        data_list =[sample]
                        results_for_file=bySampleResults[sample]
                        if (len(results_for_file.keys()) < 1):
                                data_list.append("Analysis fail. No alleles detected")
                        else:

                                for primer_set in primer_sets:
                                        for loci in panel_info[primer_set].keys():

                                                if loci in results_for_file:
                                                        allele_calls=results_for_file[loci].alleles_called
                                                else:
                                                        allele_calls=["" for x in range(0, expected_space_for_calls)]

                                                for i in range(0, expected_space_for_calls):

                                                        if (i <len(allele_calls))   :
                                                                #data_list.append(loci + ":" + str(allele_calls[i]))
                                                                data_list.append(str(allele_calls[i]))

                                                        else:
                                                                #data_list.append(loci + ":" + "-")
                                                                data_list.append("-")

                        data_line="\t".join(data_list)
                        f.write(data_line+ "\n")
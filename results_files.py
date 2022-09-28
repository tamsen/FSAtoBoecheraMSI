import os
from datetime import datetime

def WriteResults(outputDir, data):

        resultsFile = os.path.join(outputDir,"Results.csv")
        now = datetime.now()
        day = now.strftime("%d/%m/%Y")
        time = now.strftime("%H:%M:%S")
        data_string = ",".join(data)
        time_stamp_string = ",".join([day,time])

        with open(resultsFile, 'a') as f:
           f.write(time_stamp_string + "," + data_string + "\n")

def ConsoliteByFileResultsToBySampleResults(results_by_file,panel_info ):

        bySampleResults={}
        primer_sets = panel_info.keys()

        for file in results_by_file.keys():

                base = os.path.basename(file).split(".")[0].split("_")[0]
                for primer_set in primer_sets:
                        base = base.replace(primer_set,"")

                sampleNameGuess=base

                if not(sampleNameGuess in bySampleResults):
                        bySampleResults[sampleNameGuess] = results_by_file[file]
                else:
                        for loci in results_by_file[file].keys():
                                bySampleResults[sampleNameGuess][loci] = results_by_file[file][loci]

        return bySampleResults

def WriteSummaryFile(outputDir, results_by_file, panel_info ):

        now = datetime.now()
        day = now.strftime("%d_%m_%Y")
        time = now.strftime("%H_%M_%S")
        time_stamp_string = "_".join([day,time])
        summaryFile = os.path.join(outputDir,"Summary_"+ time_stamp_string +".csv")

        primer_sets=panel_info.keys()
        expected_space_for_calls = 5

        header1_data=["PrimerSet->"]
        header2_data=["Loci->"]
        for primer_set in primer_sets:
                for loci in panel_info[primer_set].keys():
                        for i in range(0,expected_space_for_calls ):
                                header1_data.append(primer_set)
                                header2_data.append(loci)

        header1 = ",".join([str(p) for p in header1_data])
        header2 = ",".join([str(p) for p in header2_data])

        bySampleResults = ConsoliteByFileResultsToBySampleResults(results_by_file, panel_info)

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
                                                        allele_calls=results_for_file[loci]
                                                else:
                                                        allele_calls=["" for x in range(0, expected_space_for_calls)]

                                                for i in range(0, expected_space_for_calls):

                                                        if (i <len(allele_calls))   :
                                                                #data_list.append(loci + ":" + str(allele_calls[i]))
                                                                data_list.append(str(allele_calls[i]))

                                                        else:
                                                                #data_list.append(loci + ":" + "-")
                                                                data_list.append("-")

                        data_line=",".join(data_list)
                        f.write(data_line+ "\n")
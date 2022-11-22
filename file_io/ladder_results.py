import os

from signal_processing.elastic_ladder_analysis import GLOBAL_Liz500


def write_ladder_peaks(outputDir, fsa_results_by_file_name):

    ladder_file_name = os.path.join(outputDir, "LadderPeaks.csv")
    peak_names = [str(n) for n in GLOBAL_Liz500]
    header_x = ["peak x's:"] + peak_names
    header_y = ["peak y's:"] + peak_names

    with open(ladder_file_name, 'a') as f:

        f.write("\n")
        f.write("\t".join(header_x) + "\n")

        for fsa_file in fsa_results_by_file_name:
            base = os.path.basename(fsa_file)
            results = fsa_results_by_file_name[fsa_file]
            data = [base]
            if results:
                ladder_peaks = results.ladder
                for peak in ladder_peaks:
                    data.append(str(peak[0]))

            f.write("\t".join(data) + "\n")

        f.write("\n")
        f.write("\t".join(header_y) + "\n")

        for fsa_file in fsa_results_by_file_name:
            base = os.path.basename(fsa_file)
            results = fsa_results_by_file_name[fsa_file]
            data = [base]
            if results:
                ladder_peaks = results.ladder
                for peak in ladder_peaks:
                    data.append(str(peak[1]))

            f.write("\t".join(data) + "\n")

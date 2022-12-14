import os
import statistics
from datetime import datetime


def allele_accuracy(called_alleles, expected_alleles):

    if len(expected_alleles) == 0:
        return -1

    #not fair to quibble about this.. there was no call in the original true data
    if expected_alleles == [-1] :
        if len(called_alleles) in [0, 1]:
            return 100.0

    #not fair to quibble about this.. there was no call in the original true data
    if expected_alleles == [0] :
        if len(called_alleles) == 0:
            return 100.0

    if len(called_alleles) == 0:
        return 0

    # dont put a null or missing allele in this calculation

    have_a_missing_allele_in_expectations = False
    null_or_missing_allele_codes = [0, -1]
    cleaned_expected_alleles = []
    for a in expected_alleles:
        if a not in null_or_missing_allele_codes:
            cleaned_expected_alleles.append(a)

    if len(cleaned_expected_alleles) == 0:
        return 100.0

    for null in null_or_missing_allele_codes:
        if null in expected_alleles:
            have_a_missing_allele_in_expectations = True

    num_expected = float(len(expected_alleles))
    num_observed = float(len(called_alleles))

    if have_a_missing_allele_in_expectations:
        diff_in_num_alleles_called = 0
    else:
        diff_in_num_alleles_called = 100.0 * (abs(num_expected - num_observed)) / max(num_expected, num_observed)

    if diff_in_num_alleles_called > 0:
        return round(100.0 - diff_in_num_alleles_called, 1)

    percent_differences = 0
    for a_exp in cleaned_expected_alleles:
        # get distance to closest called allele
        diffs = [abs(a_exp - a_obs) for a_obs in called_alleles]
        closest_it_comes = min(diffs)
        as_a_percent = 100.0 * closest_it_comes / a_exp
        percent_differences = percent_differences + as_a_percent

    return round(100.0 - percent_differences, 1)


def write_accuracy_files(outputDir, bySampleResults, panel_info):
    avg_loci_accuracy, avg_sample_accuracy = write_big_accuracy_file(outputDir, bySampleResults, panel_info)
    write_loci_accuracy_file(outputDir, avg_loci_accuracy)
    write_sample_accuracy_file(outputDir, avg_sample_accuracy)

    return avg_loci_accuracy, avg_sample_accuracy


def write_loci_accuracy_file(outputDir, avg_accuracy_for_loci_list):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])
    outFile = os.path.join(outputDir, "LociAccuracy_" + time_stamp_string + ".tsv")

    ordered_loci_list = ["ICE3", "BF20", "A1",
                         "BF11", "ICE14", "C8",
                         "BF9", "BF18", "E9",
                         "BF3", "BF19", "B6",
                         "BF15", "Bdru266", "A3"]

    header1 = "Loci\tAccuracy"

    with open(outFile, 'w') as f:
        f.write(header1 + "\n")

        for i in range(0, len(ordered_loci_list)):
            loci = ordered_loci_list[i]
            data_list = [loci, avg_accuracy_for_loci_list[i]]
            data_line = "\t".join(data_list) + "\n"
            f.write(data_line)


def write_sample_accuracy_file(outputDir, avg_sample_accuracy):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])
    outFile = os.path.join(outputDir, "SampleAccuracy_" + time_stamp_string + ".tsv")

    sample_list = avg_sample_accuracy.keys()

    header1 = "Sample\tAccuracy"

    with open(outFile, 'w') as f:
        f.write(header1 + "\n")

        for sample in sample_list:
            data_list = [sample, str(avg_sample_accuracy[sample])]
            data_line = "\t".join(data_list) + "\n"
            f.write(data_line)


def write_big_accuracy_file(outputDir, avg_accuracy_for_loci_list):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])
    outFile = os.path.join(outputDir, "LociAccuracy_" + time_stamp_string + ".tsv")

    ordered_loci_list = ["ICE3", "BF20", "A1",
                         "BF11", "ICE14", "C8",
                         "BF9", "BF18", "E9",
                         "BF3", "BF19", "B6",
                         "BF15", "Bdru266", "A3"]

    header1 = "Loci\tAccuracy"

    with open(outFile, 'w') as f:
        f.write(header1 + "\n")

        for i in range(0, len(ordered_loci_list)):
            loci = ordered_loci_list[i]
            data_list = [loci, avg_accuracy_for_loci_list[i]]
            data_line = "\t".join(data_list) + "\n"
            f.write(data_line)


def write_big_accuracy_file(outputDir, bySampleResults, panel_info):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])
    summaryFile = os.path.join(outputDir, "BigAccuracyFile_" + time_stamp_string + ".tsv")

    ordered_loci_list = ["ICE3", "BF20", "A1",
                         "BF11", "ICE14", "C8",
                         "BF9", "BF18", "E9",
                         "BF3", "BF19", "B6",
                         "BF15", "Bdru266", "A3"]

    primer_sets = panel_info.keys()
    samples = list(bySampleResults.keys())

    accuracy_by_sample = dict(zip(samples, [[] for s in samples]))
    accuracy_by_loci = dict(zip(ordered_loci_list, [[] for l in ordered_loci_list]))

    header1_data = ["PrimerSet->"]
    header2_data = ["Loci->"]
    for primer_set in primer_sets:
        for loci in panel_info[primer_set].keys():
            header1_data = header1_data + [primer_set, primer_set, primer_set]
            header2_data.append("observed " + loci)
            header2_data.append("expected " + loci)
            header2_data.append("accuracy score for " + loci)

    header2_data.append("overall sample accuracy")
    header1 = "\t".join([str(p) for p in header1_data])
    header2 = "\t".join([str(p) for p in header2_data])

    with open(summaryFile, 'w') as f:

        f.write(header1 + "\n")
        f.write(header2 + "\n")

        # for file in results_by_file:
        for sample_name in samples:

            data_list = [sample_name]
            sample_result = bySampleResults[sample_name]
            accuracies_for_this_sample = []

            for primer_set in primer_sets:
                for loci in panel_info[primer_set].keys():

                    if loci in sample_result:
                        expected_alleles = sample_result[loci].truth_data
                        called_alleles = sample_result[loci].final_alleles_called
                        accuracy_score = sample_result[loci].final_accuracy
                        accuracy_by_loci[loci].append(accuracy_score)
                        accuracies_for_this_sample.append(accuracy_score)

                        if len(called_alleles) == 0:
                            called_alleles = ["-"]

                        data_list.append(str(called_alleles))
                        data_list.append(str(expected_alleles))
                        data_list.append(str(accuracy_score))

            if len(accuracies_for_this_sample) > 0:
                sample_accuracy = statistics.mean(accuracies_for_this_sample)
            else:
                sample_accuracy = 0

            data_list.append(str(sample_accuracy))
            accuracy_by_sample[sample_name] = sample_accuracy

            data_line = "\t".join(data_list) + "\n"
            f.writelines([data_line])

        f.write("\n")
        f.write("Loci->\t" + "\t".join(ordered_loci_list) + "\n")

        avg_accuracy_for_loci_list = []
        for loci in ordered_loci_list:
            if loci in accuracy_by_loci:
                if len(accuracy_by_loci[loci]) > 0:
                    avg_accuracy_for_loci_list.append(str(statistics.mean(accuracy_by_loci[loci])))
                else:
                    avg_accuracy_for_loci_list.append("-")
            else:
                avg_accuracy_for_loci_list.append("-1")

        loci_accuracy_line = "\t".join(avg_accuracy_for_loci_list) + "\n"
        f.writelines(["Avg accuracy for loci\t" + loci_accuracy_line])
        return avg_accuracy_for_loci_list, accuracy_by_sample


def find_truth_for_this_sample(sample, truth_info):
    truth_for_this_sample = False

    if sample in truth_info.keys():
        truth_for_this_sample = truth_info[sample]
        return truth_for_this_sample

    second_try = sample.split("c")[0]
    if second_try in truth_info.keys():
        truth_for_this_sample = truth_info[second_try]
        return truth_for_this_sample

    third_try = sample.split("C")[0]
    if third_try in truth_info.keys():
        truth_for_this_sample = truth_info[third_try]

    return truth_for_this_sample

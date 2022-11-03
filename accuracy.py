import os
from datetime import datetime


def allele_accuracy(called_alleles, expected_alleles):

    # dont put a null or missing allele in this calculation

    have_a_missing_allele_in_expectations=False
    null_or_missing_allele_codes = [0, -1]
    for a in null_or_missing_allele_codes:
        if a in expected_alleles:
            expected_alleles.remove(a)
            have_a_missing_allele_in_expectations = True

    num_expected = float(len(expected_alleles))
    num_observed = float(len(called_alleles))

    if have_a_missing_allele_in_expectations:
        diff_in_num_alleles_called = 0
    else:
        diff_in_num_alleles_called = 100.0 * (abs(num_expected - num_observed)) / max(num_expected,num_observed)

    if diff_in_num_alleles_called > 0:
        return round(100.0 - diff_in_num_alleles_called,1)

    percent_differences = 0

    for a_exp in expected_alleles:
        # get distance to closest called allele
        diffs = [abs(a_exp - a_obs) for a_obs in called_alleles]
        closest_it_comes = min(diffs)
        as_a_percent = 100.0 * closest_it_comes / a_exp
        percent_differences = percent_differences + as_a_percent

    return round(100.0 - percent_differences, 1)

def assess_accuracy(outputDir, bySampleResults, panel_info):
    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])
    summaryFile = os.path.join(outputDir, "Accuracy_" + time_stamp_string + ".tsv")

    primer_sets = panel_info.keys()
    # expected_space_for_calls = 5

    header1_data = ["PrimerSet->"]
    header2_data = ["Loci->"]
    for primer_set in primer_sets:
        for loci in panel_info[primer_set].keys():
            header1_data = header1_data + [primer_set, primer_set, primer_set]
            header2_data.append("observed " + loci)
            header2_data.append("expected " + loci)
            header2_data.append("accuracy score for " + loci)

    header1 = "\t".join([str(p) for p in header1_data])
    header2 = "\t".join([str(p) for p in header2_data])

    with open(summaryFile, 'w') as f:

        f.write(header1 + "\n")
        f.write(header2 + "\n")

        # for file in results_by_file:
        for sample_name in bySampleResults.keys():

            data_list = [sample_name]
            sample_result = bySampleResults[sample_name]

            # if len(results_for_file.keys()) < 1:
            #    data_list.append("Analysis fail. No alleles detected")

            for primer_set in primer_sets:
                for loci in panel_info[primer_set].keys():

                    if loci in sample_result:
                        expected_alleles = sample_result[loci].truth_data
                        called_alleles = sample_result[loci].alleles_called
                        accuracy_score = "foo"

                        if len(called_alleles) == 0:
                            called_alleles = ["-"]

                        data_list.append(str(called_alleles))
                        data_list.append(str(expected_alleles))
                        data_list.append(str(accuracy_score))

            data_line = "\t".join(data_list) + "\n"
            f.writelines([data_line])


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

import os
from datetime import datetime


def write_query_file(outputDir, bySampleResults, print_final_calls):

    now = datetime.now()
    day = now.strftime("%d_%m_%Y")
    time = now.strftime("%H_%M_%S")
    time_stamp_string = "_".join([day, time])

    ordered_loci_list = ["ICE3", "A1", "BF20",
                         "BF11", "C8", "ICE14",
                         "BF9", "E9", "BF18",
                         "BF3", "B6", "BF19",
                         "BF15", "A3", "Bdru266"]

    how_query_likes_it = {}

    how_query_likes_it["ICE3"] = ["I3", 6]
    how_query_likes_it["BF20"] = ["B20", 4]
    how_query_likes_it["A1"] = ["A1", 3]

    how_query_likes_it["BF11"] = ["B11", 4]
    how_query_likes_it["ICE14"] = ["I14", 4]
    how_query_likes_it["C8"] = ["C8", 4]

    how_query_likes_it["BF9"] = ["B9", 4]
    how_query_likes_it["BF18"] = ["B18", 4]
    how_query_likes_it["E9"] = ["E9", 4]

    how_query_likes_it["BF3"] = ["BF3", 5]
    how_query_likes_it["BF19"] = ["BF19", 6]
    how_query_likes_it["B6"] = ["B6", 7]

    how_query_likes_it["BF15"] = ["BF15", 4]
    how_query_likes_it["Bdru266"] = ["B266", 4]
    how_query_likes_it["A3"] = ["A3", 6]

    if print_final_calls:
        summary_file = os.path.join(outputDir, "BatchQuery_FinalCallsBySample" + time_stamp_string + ".tsv")
    else:
        summary_file = os.path.join(outputDir, "BatchQuery_RawCallsBySample" + time_stamp_string + ".tsv")

    header1_data = ["Query_ID"]
    for loci in ordered_loci_list:
        instructions = how_query_likes_it[loci]
        query_loci_name = instructions[0]
        query_loci_count = instructions[1]

        for i in range(0, query_loci_count):
            header1_data.append(query_loci_name)

    header1 = "\t".join([str(p) for p in header1_data])

    with open(summary_file, 'w') as f:

        f.write(header1 + "\n")

        # for file in results_by_file:
        for sample in bySampleResults.keys():
            write_out_sample_data(bySampleResults,
                                  f, print_final_calls, ordered_loci_list, how_query_likes_it,
                                  sample)


def write_out_sample_data(bySampleResults, f, final_calls,
                          ordered_loci_list, how_query_likes_it, sample):
    data_list = [sample]
    results_for_file = bySampleResults[sample]

    if (len(results_for_file.keys()) < 1):
        data_list.append("Analysis fail. No alleles detected")
    else:

        for loci in ordered_loci_list:
            instructions = how_query_likes_it[loci]
            write_out_loci_data(data_list, instructions, final_calls, loci, results_for_file)

    data_line = "\t".join(data_list)
    f.write(data_line + "\n")


def write_out_loci_data(data_list, instructions, final_calls, loci, results_for_file):
    expected_space_for_calls = instructions[1]

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
            data_list.append("")

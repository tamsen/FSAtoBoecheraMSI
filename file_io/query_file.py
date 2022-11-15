import os
from datetime import datetime
import lxml
import requests
import pandas as pd


# to download query
# http://sites.biology.duke.edu/windhamlab/files/TD21RP21_SearchResults_ASscore.xls
import requests


def post_batch_file_and_get_response(output_dir, batch_file, bySampleResults):

    URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Search_database_batch.py"
    batch_file_result = submit_file_query(batch_file, URL1)
    destination = os.path.join(output_dir, "BatchQuery_Ouput.html")

    with open(destination, 'wb') as f:
        f.write(batch_file_result)

    results_by_sample = get_data_from_Search_database_batch_response(batch_file_result)
    # URL2 = "https://sites.biology.duke.edu/windhamlab/files/" + sample + "_SearchResults_ASscore.xls"

    # to download query
    for sample in bySampleResults:

        results_by_loci = bySampleResults[sample]

        for loci in results_by_loci:
            species_determination_data = results_by_sample[sample]
            results_by_loci[loci].set_BMW_determination(species_determination_data)

def get_data_from_Search_database_batch_response(response):

        tables = pd.read_html(response)  # Returns list of all tables on page
        results_by_sample_name={}
        num_samples=len(tables)

        for i in range(0,num_samples):

            table=tables[i]
            DNA_extractions=table['DNA extraction no.']
            species = table['Species name']
            similarity_score  = table['Allele similarity score']
            query_sample_name = DNA_extractions[0].split(" ")[1].strip()  #ie, "query: TD21RP25	"
            closest_sample_name= DNA_extractions[1]
            closest_species = species[1]
            similarity_score = similarity_score[1]
            results_by_sample_name[query_sample_name] = [closest_sample_name,closest_species,similarity_score]

        return results_by_sample_name

def submit_plain_query(URL):
    response = requests.post(URL)
    return response.content

def submit_file_query(query_filename, URL):
    with open(query_filename, 'rb') as f:
        files = {'file': f}
        print("Posting ", query_filename)
        response = requests.post(URL, files=files)

    return response.content


def write_query_file(outputDir, bySampleResults, print_final_calls): #print final calls or high-std calls
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
        batch_file = os.path.join(outputDir, "BatchQuery_FinalCallsBySample" + time_stamp_string + ".tsv")
    else:
        batch_file = os.path.join(outputDir, "BatchQuery_HighStandardCallsBySample" + time_stamp_string + ".tsv")

    header1_data = ["Query_ID"]
    for loci in ordered_loci_list:
        instructions = how_query_likes_it[loci]
        query_loci_name = instructions[0]
        query_loci_count = instructions[1]

        for i in range(0, query_loci_count):
            header1_data.append(query_loci_name)

    header1 = "\t".join([str(p) for p in header1_data])

    with open(batch_file, 'w') as f:

        f.write(header1 + "\n")

        # for file in results_by_file:
        for sample in bySampleResults.keys():
            write_out_sample_data(bySampleResults,
                                  f, print_final_calls, ordered_loci_list, how_query_likes_it,
                                  sample)

    return batch_file


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
            allele_calls = results_for_file[loci].get_high_standard_allele_calls()
    else:
        allele_calls = ["" for x in range(0, expected_space_for_calls)]
    for i in range(0, expected_space_for_calls):

        if (i < len(allele_calls)):
            data_list.append(str(allele_calls[i]))

        else:
            data_list.append("")

import os.path
from file_io import text_file_readers, xml_file_readers, results_files, query_file, ladder_results
import accuracy
from visualization import per_sample_visuals
import fsa_file_processor
import log


def process_directory(version_info, all_results_by_file, output_folder_inside_data_folder, panel_info, path,
                      results_specific_to_this_subfolder, truth_info, rules):
    for file in os.listdir(path):
        if file.endswith(".fsa"):
            fsa_file = os.path.join(path, file)

            FSA_file_results = fsa_file_processor.process_fsa_file(fsa_file,panel_info, rules,
                                                output_folder_inside_data_folder)

            if FSA_file_results:
                all_results_by_file[fsa_file] = FSA_file_results.MSI_loci_results_by_loci
                results_specific_to_this_subfolder[fsa_file] = FSA_file_results
            else:
                all_results_by_file[fsa_file] = False
                results_specific_to_this_subfolder[fsa_file] = False

    ladder_results.write_ladder_peaks(output_folder_inside_data_folder,results_specific_to_this_subfolder)

    #by_sample_results [sample_name][loci] = msi_results_for_loci
    by_sample_results = results_files.consolidate_by_file_results_to_by_sample_results(
        results_specific_to_this_subfolder, panel_info, truth_info)

    log.write_to_log("writing summary results for folder" + output_folder_inside_data_folder)

    results_files.write_summary_file(output_folder_inside_data_folder,
                                     by_sample_results, panel_info, True)

    results_files.write_summary_file(output_folder_inside_data_folder,
                                     by_sample_results, panel_info, False)

    avg_loci_accuracy, avg_sample_accuracy = accuracy.write_accuracy_files(output_folder_inside_data_folder,
                                                                           by_sample_results, panel_info)

    final_calls_batch_file = query_file.write_query_file(output_folder_inside_data_folder,
                                             by_sample_results, False, False)

    high_standard_batch_file = query_file.write_query_file(output_folder_inside_data_folder,
                                             by_sample_results, True, False)

    query_file.post_batch_file_and_get_response(output_folder_inside_data_folder, high_standard_batch_file,
                                                by_sample_results)

    species_determination_batch_file = query_file.write_query_file(output_folder_inside_data_folder,
                                             by_sample_results, True, True)

    per_sample_visuals.write_per_sample_summary_plots(version_info, output_folder_inside_data_folder,
                                                      by_sample_results)

    return by_sample_results, avg_loci_accuracy, avg_sample_accuracy

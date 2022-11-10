import os.path
from file_io import text_file_readers, xml_file_readers, results_files
import accuracy
from visualization import per_sample_visuals
import fsa_file_processor
import log


def process_directory(all_results_by_file, output_folder_inside_data_folder, panel_info, path,
                      results_specific_to_this_subfolder, truth_info):
    for file in os.listdir(path):
        if file.endswith(".fsa"):
            fsa_file = os.path.join(path, file)

            FSA_file_results = fsa_file_processor.process_fsa_file(fsa_file,
                                                                   panel_info, output_folder_inside_data_folder)

            all_results_by_file[fsa_file] = FSA_file_results.MSI_loci_results_by_loci
            results_specific_to_this_subfolder[fsa_file] = FSA_file_results

    by_sample_results = results_files.consolidate_by_file_results_to_by_sample_results(
        results_specific_to_this_subfolder, panel_info, truth_info)

    results_files.write_summary_file(output_folder_inside_data_folder,
                                     by_sample_results, panel_info)

    log.write_to_log("writing summary results for folder" + output_folder_inside_data_folder)

    per_sample_visuals.write_per_sample_summary_plots(output_folder_inside_data_folder,
                                                      by_sample_results)

    avg_loci_accuracy, avg_sample_accuracy = accuracy.write_accuracy_files(output_folder_inside_data_folder,
                                                                           by_sample_results, panel_info)

    return by_sample_results, avg_loci_accuracy, avg_sample_accuracy

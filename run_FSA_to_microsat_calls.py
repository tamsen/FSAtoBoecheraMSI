import os.path
import sys

from file_io import text_file_readers, xml_file_readers, results_files
import accuracy
import per_sample_visuals
import fsa_file_processor
import log


#code snippets from
#https://github.com/trmznt/fatools/

# Convention, https://peps.python.org/pep-0008/#function-and-variable-names
#Function names should be lowercase, with words separated by underscores as necessary to improve readability.

def greet():
    print('Tool to convert ABI FSA trace files to allele calls')

def usage():
    print('Usage:')
    print('\t%s command [options]' % sys.argv[0])
    sys.exit(0)


def main():

    greet()

    FSA_File_list = sys.argv[1]
    Panel_File = sys.argv[2]
    truth_info = {}

    if (len(sys.argv)>3):
        truth_file = sys.argv[3]
        truth_info = xml_file_readers.read_truth_data(truth_file)

    output_dir="./tmp/"
    if not(os.path.exists(output_dir)):
            os.makedirs(output_dir)

    log.write_start_to_log(output_dir)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)

    paths_to_process = text_file_readers.readInputFile(FSA_File_list)
    panel_info = xml_file_readers.readPanelXml(Panel_File)

    all_results_by_file={}

    for path in paths_to_process:

        output_folder_inside_data_folder = os.path.join(path, "FSA_to_microsat_script_results")
        results_specific_to_this_subfolder={}

        if os.path.isdir(path):

            for file in os.listdir(path):
                if file.endswith(".fsa"):
                    fsa_file=os.path.join(path, file)

                    FSA_file_results= fsa_file_processor.process_fsa_file(fsa_file,
                                                                              panel_info, output_folder_inside_data_folder)

                    all_results_by_file[fsa_file] = FSA_file_results.MSI_loci_results_by_loci
                    results_specific_to_this_subfolder[fsa_file] =  FSA_file_results

            bySampleResults = results_files.consolidate_by_file_results_to_by_sample_results(
                results_specific_to_this_subfolder, panel_info, truth_info)


            results_files.write_summary_file(output_folder_inside_data_folder,
                                             bySampleResults, panel_info)

            per_sample_visuals.write_per_sample_summary_plots(output_folder_inside_data_folder,
                                                            bySampleResults)

            accuracy.assess_accuracy(output_folder_inside_data_folder,
                                            bySampleResults, panel_info)

        else:
            final_calls_by_loci = fsa_file_processor.process_fsa_file(path, panel_info, output_dir)
            all_results_by_file[path] = final_calls_by_loci

    # will write data to the "tmp" folder or what ever is given as the output folder
    # results_files.write_summary_file(output_dir, results_by_file, panel_info)
    # accuracy_assessment.assess_accuracy(output_dir, results_by_file, panel_info, truth_info)

    log.write_end_to_log()

main()

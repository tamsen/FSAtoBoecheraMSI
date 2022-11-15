import os.path
import version
import sys
from datetime import datetime

import fsa_directory_processor
from file_io import text_file_readers, xml_file_readers, results_files
import accuracy
from visualization import per_sample_visuals
import fsa_file_processor
import log


# code snippets from
# https://github.com/trmznt/fatools/

# Convention, https://peps.python.org/pep-0008/#function-and-variable-names
# Function names should be lowercase, with words separated by underscores as necessary to improve readability.

def greet(version_info):
    print(version_info.version_num)
    print('Tool to convert ABI FSA trace files to allele calls')
    print('Example input: ' + "./data/FSAlist.txt ./data/Panel.xml ./data/TruthData.xml")

def usage():
    print('Usage:')
    print('\t%s command [options]' % sys.argv[0])
    sys.exit(0)


def main():

    version_info = version.version_info()
    greet(version_info)

    output_dir = sys.argv[1]
    FSA_File_list = sys.argv[2]
    Panel_File = sys.argv[3]
    truth_info = {}

    if (len(sys.argv) > 4):
        truth_file = sys.argv[4]
        truth_info = xml_file_readers.read_truth_data(truth_file)

    log.write_start_to_log(output_dir, version_info)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)

    paths_to_process = text_file_readers.readInputFile(FSA_File_list)
    panel_info = xml_file_readers.readPanelXml(Panel_File)
    all_results_by_file = {}

    for path in paths_to_process:

        now = datetime.now()
        day = now.strftime("%d_%m_%Y")
        time = now.strftime("%H_%M_%S")
        time_stamp_string = "_".join([day, time])

        output_folder_inside_data_folder = os.path.join(path, "FSA_to_microsat_script_results_" +
                                                        version_info.version_num +
                                                        time_stamp_string )
        results_specific_to_this_subfolder = {}

        if os.path.isdir(path):

            fsa_directory_processor.process_directory(version_info, all_results_by_file, output_folder_inside_data_folder, panel_info, path,
                              results_specific_to_this_subfolder, truth_info)

        else:
            print("Please use directories, not individual FSA files.")

    log.write_end_to_log()


main()

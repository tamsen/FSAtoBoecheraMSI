import os.path
import version
import sys
import input_parser
from datetime import datetime
import fsa_directory_processor
from file_io import text_file_readers, xml_file_readers, results_files

import log

#o=/home/tamsen/Data/Eton i=./data/FSAlist.txt panel=./data/Panel_mw.xml truth=./data/TruthData.xml rules="MW"
#o=/home/tamsen/Data/Eton i=./data/FSAlist.txt panel=./data/Panel_td.xml truth=./data/TruthData.xml rules="TD"

# code snippets from
# https://github.com/trmznt/fatools/

# Convention, https://peps.python.org/pep-0008/#function-and-variable-names
# Function names should be lowercase, with words separated by underscores as necessary to improve readability.

def greet(version_info):
    print(version_info.version_num)
    print('Tool to convert ABI FSA trace files to allele calls')
    print('Example input: ' )
    input_parser.show_example_usage()

def usage():
    print('Usage:')
    input_parser.show_example_usage()
    sys.exit(0)


def main():

    version_info = version.version_info()
    greet(version_info)
    [output_dir, FSA_File_list, Panel_File, truth_file, rules, ladder] = input_parser.do_parsing(sys.argv)
    file_path = os.path.realpath(__file__)
    internal_data_path=os.path.join(os.path.dirname(file_path),'data')

    if truth_file:
        truth_info = xml_file_readers.read_truth_data(truth_file)
    else:
        truth_info = {}

    log.write_start_to_log(output_dir, version_info)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)
    log.write_to_log('Using post-processing rule: ' + rules + "'s rules.")
    paths_to_process = text_file_readers.readInputFile(FSA_File_list)
    panel_info = xml_file_readers.readPanelXml(Panel_File)
    ladder_info = xml_file_readers.readLadderXml(os.path.join(internal_data_path,"Ladders.xml"))
    all_results_by_file = {}

    for path in paths_to_process:

        now = datetime.now()
        day = now.strftime("%d_%m_%Y")
        time = now.strftime("%H_%M_%S")
        time_stamp_string = "_".join([day, time])

        eton_order_num_dir=path.split("/")[-1]
        eton_order_num_dir_splat=eton_order_num_dir.split("_")
        if len(eton_order_num_dir_splat) > 1:
            eton_order_num = eton_order_num_dir_splat[-1] + "_"
        else:
            eton_order_num = ""

        output_folder_inside_data_folder = os.path.join(path, eton_order_num + "FSA_to_microsat_script_results__" +
                                                        version_info.version_num + "__" +
                                                        time_stamp_string )
        results_specific_to_this_subfolder = {}

        if os.path.isdir(path):

            fsa_directory_processor.process_directory(version_info, all_results_by_file, output_folder_inside_data_folder, panel_info, path,
                              results_specific_to_this_subfolder, truth_info,
                              rules, ladder_info[ladder], ladder)

        else:
            print("Please use directories, not individual FSA files.")

    log.write_end_to_log()


main()

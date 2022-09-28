import os.path
import sys
import input_file_readers
import fsa_file_processor
import log
import results_files

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
    output_dir="./tmp/"
    log.write_start_to_log(output_dir)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)

    paths_to_process = input_file_readers.readInputFile(FSA_File_list)
    panel_info = input_file_readers.readPanelXml(Panel_File)
    results_by_file={}

    for path in paths_to_process:

        if os.path.isdir(path):

            for file in os.listdir(path):
                if file.endswith(".fsa"):
                    fsa_file=os.path.join(path, file)
                    final_calls_by_loci = fsa_file_processor.processFSAfile(fsa_file, panel_info, output_dir)
                    results_by_file[fsa_file] = final_calls_by_loci
        else:
            final_calls_by_loci = fsa_file_processor.processFSAfile(path, panel_info, output_dir)
            results_by_file[path] = final_calls_by_loci

    results_files.WriteSummaryFile(output_dir, results_by_file, panel_info)

    log.write_end_to_log()

main()

import sys
import InputFileReaders
import tracefileprocessor

#code snippets from
#https://github.com/trmznt/fatools/

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

    print('Command Arguments Given: %s' % sys.argv)

    files_to_process = InputFileReaders.readInputFile(FSA_File_list)
    panel_info = InputFileReaders.readPanelXml(Panel_File)

    for file in files_to_process:
        tracefileprocessor.processFSAfile(file, panel_info)

main()

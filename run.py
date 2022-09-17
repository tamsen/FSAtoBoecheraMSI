import sys, os
import argparse
import filereaders
import tracefileprocessor

#code snippets from
#https://github.com/trmznt/fatools/

def greet():
    print('Tool to convert ABI FSA trace files to allele calls')

def usage():
    print('Usage:')
    print('\t%s command [options]' % sys.argv[0])
    sys.exit(0)

def init_argparser():

    parser = argparse.ArgumentParser( 'convert' )
    parser.add_argument('infiles', nargs='+')

    return parser

def main():

    greet()

    command = sys.argv[1]
    opt_args = sys.argv[2:]

    print('Running command: %s' % command)

    #try:

        #set up imports..

    #except ImportError:
    #    print('Cannot import script name: %s' % command)
    #    raise

    parser = init_argparser()
    args = parser.parse_args(opt_args)
    #M.main( args )

    files_to_process = filereaders.readInputFile("./data/FSAlist.txt")
    panel_info = filereaders.readPanelFile("./data/Panel.txt")

    for file in files_to_process:
        tracefileprocessor.processFSAfile(file, panel_info)

main()

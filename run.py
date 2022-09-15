import sys, os
import argparse
import FSAreader
import analysis
import visuals

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

    print("in main")

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

    trace_data, all_collected_data = FSAreader.readFSAFile("./data/TD21DG24PS1c10_C07.fsa")

    # get peaks of ladder
    sixteen_peaks = analysis.getLadderPeaks(all_collected_data)
    #visuals.plotLadder(all_collected_data, smoothed_ladder_data, indices, peak_heights)
    visuals.plotRawTraceByColor(all_collected_data)



main()

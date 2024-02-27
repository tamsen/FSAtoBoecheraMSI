# FSA to Boechera MSI


# Quick Start

Command syntax:  python run_FSA_to_microsat_calls.py Files.txt Panel.xml

"Files.txt" is a txt file with th list of FSA files you want to process.
And "Panel.xml" decribes your multiplex primer set.
(Sorry, DNA ladder currenly hard coded to LIZ 500, thats all we are using right now.)

Example, using the sample input files included in the "data" folder in the repo

tamsen@tamsen:~/Git/FSAtoBoecheraMSI/FSAtoBoecheraMSI$ /venv/bin/python run_FSA_to_microsat_calls.py data/FSAlist.txt data/Panel.xml

~ Running Tamsen's raw data ~
run_FSA_to_microsat_calls.py o=path_to_output_folder i=/path_to/FSAlist.txt panel=/path_to/Panel_td.xml truth=/path_to/TruthData.xml rules="TD" ladder="Liz500"

# Running Michael W's raw data
run_FSA_to_microsat_calls.py o=path_to_output_folder i=/path_to/FSAlist.txt panel=/path_to/Panel_mw.xml truth=/path_to/TruthData.xml rules="MW" ladder="Liz500"

# Running Teri's raw data
run_FSA_to_microsat_calls.py o=path_to_output_folder i=/path_to/FSAlist.txt panel=/path_to/Panel_mw.xml truth=/path_to/TruthData.xml rules="TM" ladder="400HD"


# Code overview

The code loops through each file:

For each file, it reads in the ladder channel and data channels in the FSA file.
  The  TraceAnalysis.getLadderPeaks picks off where the 16 ladder peaks are in the ladder trace.
  The code builds a cubic spline interpolation mapping between the peaks in base-pair-space and the x-axis of the trace (ie, the raw distance traveled by       the fragment in the gel).

For each channel (dye color) there is a data channel trace. The code loops through each relevant data channel.
    The mapping (derived from the ladder) is used to map each data channels into basepair space (TraceAnalysis.RemapDataTrace) and returns the peaks in the data channels (Peaks_inside_loci) in base-pair-space.
    Guided by the panel data, only the peaks in the loci of interest according to the panel.xml file are kept.
    The peaks get filtered to remove noise (PeakAnalysis.PeaksToMsiCalls)
    And then get written out to file (ResultsFile.WriteResults(output_dir, data) )


# Output overview

Output is written to a "tmp" folder where you run the cmd.

Results per input file are written (appended) to "Results.csv". 

There is an output folder generated for every .fsa file. In that folder are plots for the ladder and the trace (in raw and bp space), the mapping to bp chart (sometimes its a nice fit, but sometimes its not monotonic, so a warning there. TODO put in a fix), and a zoomed-in graph of every allele call for every loci (so you can always see why the code made the decision it did).

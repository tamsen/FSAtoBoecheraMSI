import FSAreader
import TraceAnalysis
import os
import ResultsFile
import InputFileReaders
import PeakAnalysis


def processFSAfile(FSAfile, panel_info):

    print("****** Processing " + FSAfile + " **********")

    dye_to_channel_mapping, all_collected_data = FSAreader.readFSAFile(FSAfile)

    output_dir="./tmp/"
    run_name =all_collected_data["SpNm1"]
    run_folder = os.path.join(output_dir,run_name)
    if not(os.path.exists(run_folder )):
            os.makedirs(run_folder )


    relevant_loci = InputFileReaders.figure_out_loci_from_run_name(panel_info, run_name)

    if (relevant_loci == "FAIL" ):
        print("Uh-oh!  Can't figure out what panel to use for this FSA file!!")
        print("skipping " + run_name)
        return
    else:
        print("relevant_loci=", str(relevant_loci ))


    sixteen_peaks, threshold = TraceAnalysis.getLadderPeaks(run_folder, all_collected_data)

    mappingFxn, left_domain_limit, right_domain_limit = \
        TraceAnalysis.buildInterpolationdBasedOnLadder(run_folder, sixteen_peaks)

    TraceAnalysis.RemapLadder(run_folder, all_collected_data, mappingFxn,
                              left_domain_limit, right_domain_limit,
                              sixteen_peaks, threshold)


    # Channels we care about are ones with dyes in our panel.
    channels = set([loci_info_dict["dye"] for loci_info_dict in relevant_loci.values()])

    for channel in channels:
        Peaks_inside_loci, trace_x_new, trace_y_new  = TraceAnalysis.RemapDataTrace(run_folder,
                                                                                    relevant_loci,  #ie, the loci for this primer set
                                                                                    all_collected_data, mappingFxn,
                                                                                    left_domain_limit, right_domain_limit,
                                                                                    sixteen_peaks, dye_to_channel_mapping[channel])



        for loci in Peaks_inside_loci:
            # convert peak calls to MSI calls.
            print("loci " + loci)
            MSI_calls = PeakAnalysis.PeaksToMsiCalls(
                Peaks_inside_loci[loci],trace_x_new, trace_y_new, threshold)

            allele_calls_for_loci = [str(x[0]) for x in MSI_calls ]
            data = [ FSAfile, loci ] + allele_calls_for_loci
            ResultsFile.WriteResults(output_dir, data)

    print(run_name + " completed ")
import fsa_file_reader
import trace_analysis
import os
import results_files
import input_file_readers
import peak_analysis
import visuals
import log


def processFSAfile(FSAfile, panel_info,output_dir):

    log.write_to_log("****** Processing " + FSAfile + " **********")

    dye_to_channel_mapping, all_collected_data = fsa_file_reader.readFSAFile(FSAfile)

    run_name =all_collected_data["SpNm1"]
    run_folder = os.path.join(output_dir,run_name)
    if not(os.path.exists(run_folder )):
            os.makedirs(run_folder )


    relevant_loci = input_file_readers.figure_out_loci_from_run_name(panel_info, run_name)

    if (relevant_loci == "FAIL" ):
        log.WriteErrorToLog("Uh-oh!  Can't figure out what panel to use for this FSA file!!")
        log.WriteErrorToLog("Quitting " + run_name)
        return {}
    else:
        log.write_to_log("relevant_loci: " + str(relevant_loci))

    gotLadderPeaks = trace_analysis.getLadderPeaks(run_folder, all_collected_data)

    if not(gotLadderPeaks):
        log.write_to_log("**** Processing " + FSAfile + " failed ********")
        return {}

    else:
        sixteen_peaks, threshold =  gotLadderPeaks

    mapping_fxn, left_domain_limit, right_domain_limit = \
        trace_analysis.buildInterpolationdBasedOnLadder(run_folder, sixteen_peaks)

    trace_analysis.RemapLadder(run_folder, all_collected_data, mapping_fxn,
                               left_domain_limit, right_domain_limit,
                               sixteen_peaks, threshold)


    # Channels we care about are ones with dyes in our panel.
    channels = set([loci_info_dict["dye"] for loci_info_dict in relevant_loci.values()])
    final_calls_by_loci={}

    for channel in channels:
        peaks_inside_loci, trace_x_new, trace_y_new = trace_analysis.RemapDataTrace(run_folder,
                                                                                     relevant_loci,  #ie, the loci for this primer set
                                                                                     all_collected_data, mapping_fxn,
                                                                                     left_domain_limit, right_domain_limit,
                                                                                     sixteen_peaks, dye_to_channel_mapping[channel])


        for loci in peaks_inside_loci:

            # convert peak calls to MSI calls.
            unfiltered_peaks_in_loci=peaks_inside_loci[loci]
            log.write_to_log("Processing loci " + loci)
            log.write_to_log("Unfiltered peaks found in loci range: " + str(unfiltered_peaks_in_loci))

            MSI_calls = peak_analysis.peaks_to_msi_calls(unfiltered_peaks_in_loci, trace_x_new, trace_y_new, threshold)

            #make a zoomed-in plot JUST around the MSI call
            for final_call in MSI_calls:

                final_call_x = final_call[0]
                final_call_y = final_call[1]
                domain = [final_call_x-20, final_call_x+20]
                plot_prefix = "Loci" + loci + "CallAt" + str(final_call_x )
                visuals.plotRemappedTrace(run_folder, trace_x_new, trace_y_new,
                                          [final_call_x], [final_call_y], threshold,
                              "wavelength",str( dye_to_channel_mapping[channel]), channel,
                              plot_prefix, plot_domain=domain)

            #get the results ready to print to file
            allele_calls_for_loci = [str(x[0]) for x in MSI_calls ]
            data = [ FSAfile, loci ] + allele_calls_for_loci
            results_files.WriteResults(output_dir, data)
            final_calls_by_loci[loci] = allele_calls_for_loci
            log.write_to_log("final calls for loci " + loci + ": " + str(allele_calls_for_loci))

    log.write_to_log("**** Processing " + FSAfile + " completed  ********")

    return final_calls_by_loci
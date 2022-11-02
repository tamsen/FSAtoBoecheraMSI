import fsa_file_reader
import trace_analysis
import os
import results_files
import xml_file_readers
import peak_analysis
import visuals
import log
from results_for_fsa_file import FSA_File_Results, MSI_loci_results


def process_fsa_file(fsa_file, panel_info, output_dir):

    log.write_to_log("****** Processing " + fsa_file + " **********")
    dye_to_channel_mapping, all_collected_data = fsa_file_reader.readFSAFile(fsa_file)
    run_name =all_collected_data["SpNm1"]
    run_folder = os.path.join(output_dir,run_name)
    if not(os.path.exists(run_folder )):
            os.makedirs(run_folder )


    relevant_loci = xml_file_readers.figure_out_loci_from_run_name(panel_info, run_name)

    if (relevant_loci == False):
        log.write_error_to_log("Uh-oh!  Can't figure out what panel to use for this FSA file!!")
        log.write_error_to_log("Quitting " + run_name)
        data = [fsa_file, "panel problem", "Can't figure out what panel to use for this FSA file!!"]
        results_files.write_results(output_dir,data)
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        return {}
    else:
        log.write_to_log("relevant_loci: " + str(relevant_loci))

    gotLadderPeaks = trace_analysis.getLadderPeaks(run_folder, all_collected_data)

    if not(gotLadderPeaks):
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        data = [fsa_file, "panel problem", "Couldn't get the right number of ladder peaks!!"]
        results_files.write_results(output_dir, data)

        return {}

    else:
        sixteen_peaks, threshold, ladder_plot_data =  gotLadderPeaks

    ladder_worked = trace_analysis.build_interpolation_based_on_ladder(run_folder, sixteen_peaks)

    if not (ladder_worked):
        data = [fsa_file, "panel problem", "Ladder failed monotonicity!!"]
        results_files.write_results(output_dir, data)
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        return {}

    else:
        mapping_fxn, left_domain_limit, right_domain_limit, mapping_plot_data_spline,mapping_plot_data_linear, \
            = ladder_worked

    trace_analysis.remap_ladder(run_folder, all_collected_data, mapping_fxn, left_domain_limit, right_domain_limit,
                                sixteen_peaks, threshold)


    # Channels we care about are ones with dyes in our panel.
    channels = set([loci_info_dict["dye"] for loci_info_dict in relevant_loci.values()])
    #all_loci_plot_data=[]

    final_calls_by_loci={}


    for channel in channels:

        peaks_inside_loci, trace_x_new, trace_y_new = trace_analysis.remap_data_trace(run_folder, relevant_loci,
                                                                                      all_collected_data, mapping_fxn,
                                                                                      left_domain_limit,
                                                                                      right_domain_limit, sixteen_peaks,
                                                                                      dye_to_channel_mapping[channel])


        for loci in peaks_inside_loci:

            # convert peak calls to MSI calls.
            unfiltered_peaks_in_loci=peaks_inside_loci[loci]
            log.write_to_log("Processing loci " + loci)
            log.write_to_log("Unfiltered peaks found in loci range: " + str(unfiltered_peaks_in_loci))

            channel_specific_threshold = trace_analysis.get_threshold_for_trace(trace_y_new)
            MSI_calls_for_loci = peak_analysis.peaks_to_msi_calls(unfiltered_peaks_in_loci, trace_x_new,
                                                         trace_y_new,channel_specific_threshold )

            #make a zoomed-in plot JUST around the MSI call
            for final_call in MSI_calls_for_loci:

                final_call_x = final_call[0]
                final_call_y = final_call[1]
                domain = [final_call_x-20, final_call_x+20]
                plot_prefix = "Loci" + loci + "CallAt" + str(final_call_x )



                visuals.plot_remapped_trace(run_folder, trace_x_new, trace_y_new, [final_call_x], [final_call_y],
                                            threshold, "wavelength", str(dye_to_channel_mapping[channel]), channel,
                                            plot_prefix,plot_domain= domain)

            if (len(MSI_calls_for_loci) > 0):
                whole_loci_domain = [MSI_calls_for_loci[0][0]-20,domain[1]]

                loci_specific_plot_data=[trace_x_new, trace_y_new,
                                          [call[0] for call in MSI_calls_for_loci],
                                          [call[1] for call in MSI_calls_for_loci],
                                           loci, whole_loci_domain, channel]
            else:
                where_loci_should_be=relevant_loci[loci]["length"]  #even though we didnt find them..
                whole_loci_domain = [where_loci_should_be[0],where_loci_should_be[-1]]
                loci_specific_plot_data=[trace_x_new, trace_y_new,
                                          [call[0] for call in MSI_calls_for_loci],
                                          [call[1] for call in MSI_calls_for_loci],
                                           loci, whole_loci_domain, channel]


            #get the results ready to print to file
            allele_calls_for_loci = [str(x[0]) for x in MSI_calls_for_loci ]
            data = [fsa_file, loci] + allele_calls_for_loci
            results_files.write_results(output_dir, data)

            loci_results=MSI_loci_results(allele_calls_for_loci,
                                          loci_specific_plot_data, ladder_plot_data,
                                          [mapping_plot_data_spline, mapping_plot_data_linear])

            final_calls_by_loci[loci]=loci_results

            log.write_to_log("final calls for loci " + loci + ": " + str(allele_calls_for_loci))

    log.write_to_log("**** Processing " + fsa_file + " completed  ********")
    FSA_file_results = FSA_File_Results(final_calls_by_loci)

    return FSA_file_results
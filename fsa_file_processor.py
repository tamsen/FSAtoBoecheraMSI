import os
from file_io import xml_file_readers, results_files, fsa_file_reader
from signal_processing import peak_analysis, trace_analysis, ladder_analysis, mw_wisdom, shared
from fsa_file_results import FSA_File_Results
from loci_results import loci_results
from signal_processing.ladder_analysis import Mapping_Info
from visualization import per_file_visuals
import log


def process_fsa_file(fsa_file, panel_info, output_dir):
    log.write_to_log("****** Processing " + fsa_file + " **********")
    dye_to_channel_mapping, all_collected_data = fsa_file_reader.readFSAFile(fsa_file)
    run_name = all_collected_data["SpNm1"]
    run_folder = os.path.join(output_dir, run_name)
    if not (os.path.exists(run_folder)):
        os.makedirs(run_folder)

    relevant_loci = xml_file_readers.figure_out_loci_from_file_name(panel_info, fsa_file)

    if (relevant_loci == False):
        log.write_error_to_log("Uh-oh!  Can't figure out what panel to use for this FSA file!!")
        log.write_error_to_log("Quitting " + run_name)
        data_string = [fsa_file, "panel problem", "Can't figure out what panel to use for this FSA file!!"]
        results_files.write_results(output_dir, data_string)
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        return False
    else:
        log.write_to_log("relevant_loci: " + str(relevant_loci))

    mapping_attempt_worked = use_the_ladder_to_make_a_mapping(all_collected_data, fsa_file,
                                                              output_dir, run_folder, run_name, -1)
    retry_needed = check_if_retry_is_worth_it(mapping_attempt_worked)

    if retry_needed:
        mapping_attempt_worked = use_the_ladder_to_make_a_mapping(all_collected_data, fsa_file,
                                                             output_dir, run_folder, run_name, 50)

        retry_needed = check_if_retry_is_worth_it(mapping_attempt_worked)
        if retry_needed:
            mapping_attempt_worked = use_the_ladder_to_make_a_mapping(all_collected_data, fsa_file,
                                                             output_dir, run_folder, run_name, 25)
    if not mapping_attempt_worked:  # still!
        log.write_to_log("getting ladder peaks failed")
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        return False

    ladder_plot_data, mapping_function, sixteen_peaks, threshold = mapping_attempt_worked

    trace_analysis.remap_ladder(run_folder, all_collected_data, mapping_function, sixteen_peaks, threshold)

    # Channels we care about are ones with dyes in our panel.
    channels = set([loci_info_dict["dye"] for loci_info_dict in relevant_loci.values()])

    final_calls_by_loci = {}
    for channel in channels:

        threshold_multiplier = 0.5
        # best parameters for the trace data: (fatter, messier spikes)
        # kernel of 20, distance between peaks 20, min_peak_width 10;threshold_multiplier .5

        # Originally used a lot of smoothing, but backed off to better match michael w
        # peak_calling_parameters = shared.peak_calling_parameters(30, 20, 20, 10, threshold_multiplier)
        peak_calling_parameters = shared.peak_calling_parameters(30, 10, 3, 1, threshold_multiplier, False)

        peaks_inside_loci, trace_x_new, trace_y_new, \
        threshold_used = trace_analysis.remap_data_trace_and_call_raw_peaks(run_folder, relevant_loci,
                                                                            all_collected_data, mapping_function,
                                                                            sixteen_peaks,
                                                                            dye_to_channel_mapping[channel],
                                                                            peak_calling_parameters)

        for loci in peaks_inside_loci.keys():

            # convert peak calls to MSI calls.
            unfiltered_peaks_in_loci = peaks_inside_loci[loci]
            log.write_to_log("Processing loci " + loci)
            log.write_to_log("Unfiltered peaks found in loci range: " + str(unfiltered_peaks_in_loci))

            raw_calls = peak_analysis.peaks_to_raw_calls(unfiltered_peaks_in_loci)

            filtered_calls = peak_analysis.peaks_to_filtered_calls(raw_calls, loci)

            # rescue a loci that might be low-intensity
            # and thus falling below threshold
            if (len(raw_calls)) > 0:
                max_call_intensity = max(call[1] for call in raw_calls)
                threshold_reduction = 0.3
            else:
                max_call_intensity = 0
                threshold_reduction = 0.01

            rescue_needed = (max_call_intensity < 4 * threshold_used)
            if rescue_needed:
                rescue_parameters = shared.peak_calling_parameters(
                    100,
                    peak_calling_parameters.kernel_size,
                    peak_calling_parameters.min_distance_between_peaks,
                    peak_calling_parameters.min_distance_between_peaks,
                    threshold_multiplier * threshold_reduction, False)
                #    #False)
                #
                # threshold_multiplier * threshold_reduction)

                peaks_inside_loci, trace_x_new2, trace_y_new2, \
                threshold_used = trace_analysis.remap_data_trace_and_call_raw_peaks(run_folder, relevant_loci,
                                                                                    all_collected_data,
                                                                                    mapping_function,
                                                                                    sixteen_peaks,
                                                                                    dye_to_channel_mapping[channel],
                                                                                    rescue_parameters)

                unfiltered_peaks_in_loci = peaks_inside_loci[loci]
                raw_calls = peak_analysis.peaks_to_raw_calls(unfiltered_peaks_in_loci)

                filtered_calls = peak_analysis.peaks_to_filtered_calls(raw_calls, loci)

            final_calls = mw_wisdom.make_adjustments(filtered_calls, loci)

            # make a zoomed-in plot JUST around the MSI call
            for final_call in final_calls:
                final_call_x = final_call[0]
                final_call_y = final_call[1]
                domain = [final_call_x - 20, final_call_x + 20]
                plot_prefix = "Loci" + loci + "CallAt" + str(final_call_x)

                per_file_visuals.plot_remapped_trace(run_folder, trace_x_new, trace_y_new, [final_call_x],
                                                     [final_call_y],
                                                     threshold, "wavelength", str(dye_to_channel_mapping[channel]),
                                                     channel,
                                                     plot_prefix, plot_domain=domain)

            if (len(final_calls) > 0):
                whole_loci_domain = [final_calls[0][0] - 20, domain[1]]
            else:
                where_loci_should_be = relevant_loci[loci]["length"]  # even though we didnt find them..
                whole_loci_domain = [where_loci_should_be[0], where_loci_should_be[-1]]

            loci_specific_plot_data = [trace_x_new, trace_y_new,
                                       [call[0] for call in raw_calls],
                                       [call[1] for call in raw_calls],
                                       [call[0] for call in filtered_calls],
                                       [call[1] for call in filtered_calls],
                                       threshold_used,
                                       loci, whole_loci_domain, channel]

            # get the results ready to print to file
            raw_calls_for_loci = [x[0] for x in raw_calls]
            filtered_calls_for_loci = [x[0] for x in filtered_calls]
            allele_calls_for_loci = [x[0] for x in final_calls]
            data_string = [fsa_file, loci] + [str(x) for x in allele_calls_for_loci]
            results_files.write_results(output_dir, data_string)

            results_for_loci = loci_results(raw_calls_for_loci, filtered_calls_for_loci,
                                            allele_calls_for_loci, fsa_file,
                                            loci_specific_plot_data, ladder_plot_data,
                                            [mapping_function.mapping_plot_data_spline,
                                             mapping_function.mapping_plot_data_linear],
                                            mapping_function.ladder_state)

            final_calls_by_loci[loci] = results_for_loci

            log.write_to_log("final calls for loci " + loci + ": " + str(allele_calls_for_loci))

    log.write_to_log("**** Processing " + fsa_file + " completed  ********")
    FSA_file_results = FSA_File_Results(final_calls_by_loci,sixteen_peaks)

    return FSA_file_results


def check_if_retry_is_worth_it(mapping_attempt_worked):
    retry_needed = False
    if not mapping_attempt_worked:
        retry_needed = True
    else:
        ladder_plot_data, mapping_function, sixteen_peaks, threshold = mapping_attempt_worked
        if mapping_function.ladder_state == Mapping_Info.LadderState.Suspect:
            retry_needed = True
    return retry_needed


def use_the_ladder_to_make_a_mapping(all_collected_data, fsa_file, output_dir, run_folder, run_name,
                                     background_subtraction_window):

    gotLadderPeaks = ladder_analysis.getLadderPeaks(run_folder, run_name, all_collected_data,
                                                    background_subtraction_window)

    if not gotLadderPeaks:
        log.write_to_log("getting ladder peaks failed")
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        data_string = [fsa_file, "panel problem", "Couldn't get the right number of ladder peaks!!"]
        results_files.write_results(output_dir, data_string)
        return False

    else:
        sixteen_peaks, threshold, ladder_plot_data = gotLadderPeaks
    mapping_function = ladder_analysis.build_interpolation_based_on_ladder(run_folder, sixteen_peaks)

    if not mapping_function:
        data_string = [fsa_file, "panel problem", "Ladder failed monotonicity!!"]
        results_files.write_results(output_dir, data_string)
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        return False

    if mapping_function.ladder_state == Mapping_Info.LadderState.Suspect:
        data_string = [fsa_file, "panel problem", "Ladder too suspicious!!"]
        results_files.write_results(output_dir, data_string)
        log.write_to_log("**** Processing " + fsa_file + " failed ********")
        #return False


    return ladder_plot_data, mapping_function, sixteen_peaks, threshold

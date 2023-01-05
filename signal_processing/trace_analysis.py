from signal_processing import peak_analysis, elastic_ladder_analysis
from visualization import per_file_visuals
import log
from signal_processing.elastic_ladder import GLOBAL_Liz500

def remap_a_trace(tracedata_x_coords, tracedata_y_coords,
                  fxn, left_domain_limit, right_domain_limit):
    x_raw = []
    y_raw = []

    # get the subset of coords in the domain
    for i in range(0, len(tracedata_x_coords)):
        x = tracedata_x_coords[i]
        if x >= left_domain_limit and x <= right_domain_limit:
            x_raw.append(tracedata_x_coords[i])
            y_raw.append(tracedata_y_coords[i])

    if fxn:
        x_new = fxn(x_raw)
    else:
        x_new = x_raw

    return x_raw, x_new, y_raw  # y_raw and y_new are the same


def remap_ladder(run_folder, trace_data_dictionary, mapping_info, sixteen_peaks, threshold, ladder_channel):
    original_ladder = trace_data_dictionary[ladder_channel]
    num_data_points = len(original_ladder)
    old_x_coords = [x for x in range(0, num_data_points)]
    x_raw, x_new, y_new = remap_a_trace(old_x_coords, original_ladder, mapping_info.mapping_fxn,
                                        mapping_info.left_ladder_domain_limit, mapping_info.right_ladder_domain_limit)

    peak_xs = GLOBAL_Liz500
    peak_ys = [peak[1] for peak in sixteen_peaks]

    per_file_visuals.plot_remapped_trace(run_folder, x_new, y_new, peak_xs, peak_ys, threshold,
                                "Ladder Trace", "Ladder Trace",
                                "Ladder Trace", "Remapped")

    return True


def remap_data_trace_and_call_raw_peaks(run_folder, relevant_loci,
                                        trace_data_dictionary, mapping_info,
                                        sixteen_peaks, channel_number, peak_calling_parameters):

    channel_name = 'DATA' + str(channel_number)
    wavelength = trace_data_dictionary['DyeW' + str(channel_number)]
    channel_dye_name = trace_data_dictionary['DyeN' + str(channel_number)]
    raw_trace_data = (trace_data_dictionary[channel_name])

    highest_peaks_tup, smoothed_trace, threshold = elastic_ladder_analysis.find_top_N_Peaks(
        raw_trace_data, peak_calling_parameters, True)

    highest_peaks_tup.sort(key=lambda x: x[0])  # sort, by x's, not y's

    peak_xs = [peak[0] for peak in highest_peaks_tup]
    peak_ys = [peak[1] for peak in highest_peaks_tup]

    # plot raw trace, smoothed trace, and discovered peaks
    per_file_visuals.plotUnmappedTraceByColor(
        run_folder,
        raw_trace_data, smoothed_trace, highest_peaks_tup,
        wavelength, channel_dye_name, channel_number, "Raw")

    num_data_points = len(smoothed_trace)
    old_x_coords = [x for x in range(0, num_data_points)]
    trace_x_raw, trace_x_new, trace_y_new = remap_a_trace(old_x_coords, smoothed_trace, mapping_info.mapping_fxn,
                                                          mapping_info.left_ladder_domain_limit,
                                                          mapping_info.right_ladder_domain_limit)

    # if peaks are at the extreme end of the ladder, throw them out
    peak_x_raw, peak_x_new, peak_y_new = remap_a_trace(peak_xs, peak_ys,
                                                       mapping_info.mapping_fxn, sixteen_peaks[1][0] + 10,
                                                       sixteen_peaks[14][0] - 10)

    plot_domain = [trace_x_new[0], trace_x_new[len(trace_x_new) - 1]]
    log.write_to_log("acceptable domain based on ladder: " + str(plot_domain))
    log.write_to_log("peaks found: " + str(peak_x_new))

    per_file_visuals.plot_remapped_trace(run_folder, trace_x_new, trace_y_new,
                                         peak_x_new, peak_y_new, threshold, wavelength,
                                         channel_dye_name, channel_number, "Remapped", plot_domain)

    # Plot a zoomed-in view for each loci of interest,
    # for this dye (usually one or two loci per dye).
    # Save off the MSI calls for each plot

    Peaks_inside_loci = {}
    for loci_name in relevant_loci.keys():

        loci = relevant_loci[loci_name]
        dye_in_panel = loci["dye"]
        if dye_in_panel == channel_dye_name:
            log.write_to_log("Scanning " + dye_in_panel + " trace")

        elif (dye_in_panel == "FAM") and \
                (channel_dye_name == "6-FAM"):
            # ok, good enough.
            log.write_to_log("Scanning FAM trace")

        else:
            continue

        bp_start = loci["length"][0] - 20
        bp_end = loci["length"][1] + 20
        plot_domain = [bp_start, bp_end]

        per_file_visuals.plot_remapped_trace(run_folder, trace_x_new, trace_y_new, peak_x_new, peak_y_new, threshold, wavelength,
                                             channel_dye_name, channel_number, loci_name + "_Remapped", plot_domain)

        peaks = peak_analysis.filter_by_range(peak_x_new, peak_y_new, loci["length"])

        Peaks_inside_loci[loci_name] = peaks

    return Peaks_inside_loci, trace_x_new, trace_y_new, threshold

import numpy as np
from scipy.stats import mode
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
import per_file_visuals
import log

GLOBAL_Liz500 = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500]


def getLadderPeaks(runFolder, runName, trace_data_dictionary):
    log.write_to_log("Reading through ladder trace for " + runName)
    # 'DATA105' is the ladder channel
    ladder_trace = trace_data_dictionary['DATA105']

    # parameters_for_right_side_of_ladder
    # kernel_size=10
    # min_distance_between_peaks=5#2
    # min_peak_width=5 #2
    # threshold_multiplier=.50
    parameters_for_right_side_of_ladder = [30, 20, 50, 10, .5]
    parameters_for_left_side_of_ladder = [30, 20, 2, 1, .5]
    highest_peaks_tup, smoothed_trace, threshold = find_top_N_Peaks(
        ladder_trace, *parameters_for_right_side_of_ladder, True)

    # now, keep the best 16 starting from the right-most index.
    highest_peaks_tup.sort(key=lambda x: x[0], reverse=True)
    sixteen_peaks = highest_peaks_tup[0:16]

    # want your ladder peaks leftmost on the left! not sorted by size
    sixteen_peaks.sort(key=lambda x: x[0])

    # Minor tweak, our vendor has a spurious peak that pops up between  2000 and 2800
    # SO - if we see 4 peaks between 2000 and 2800, and we should see ony 3,
    # throw out the smallest
    sixteen_peaks = remove_known_sus_ladder_peak(sixteen_peaks, highest_peaks_tup)
    sixteen_peaks.sort(key=lambda x: x[0])

    # Another minor tweak, we need a smaller kernel to get the very first ladder point right,
    # Because it slams into over-saturation at the beginning of the gel
    best_guess_for_35_peak = fix_over_saturated_start(ladder_trace, parameters_for_left_side_of_ladder,
                                                      sixteen_peaks)
    sixteen_peaks[0]=best_guess_for_35_peak

    # bake LIZ500 into the peak tuple, for use later on.
    numLadderPeaks = len(sixteen_peaks)
    log.write_to_log("Num peaks found for ladder: " + str(numLadderPeaks))

    sixteen_peaks = [(sixteen_peaks[i][0], sixteen_peaks[i][1], GLOBAL_Liz500[i])
                     for i in range(0, numLadderPeaks)]

    ladder_plot_data = [runFolder, runName + "_LadderPlot", threshold, smoothed_trace, sixteen_peaks]
    per_file_visuals.plot_ladder(*ladder_plot_data, )

    # note, we had an index out of range here - issue with the ladder - hence the check
    if (numLadderPeaks != len(GLOBAL_Liz500)):
        log.write_to_log("There is a problem with this sample's ladder! Inspect plot. ")
        log.write_to_log("Aborting run. ")
        return False

    return sixteen_peaks, threshold, ladder_plot_data


def fix_over_saturated_start(ladder_trace, parameters_for_left_side_of_ladder, sixteen_peaks):

    #highest_peaks_tup_2, smoothed_trace_2, threshold_2 = find_top_N_Peaks(
    #    ladder_trace[0:int(len(ladder_trace) * 0.1)], *parameters_for_left_side_of_ladder, False)
    highest_peaks_tup, smoothed_trace, threshold = find_top_N_Peaks(
            ladder_trace, *parameters_for_left_side_of_ladder, False)

    original_bp35_peak = sixteen_peaks[0]
    original_bp50_peak = sixteen_peaks[1]
    furthest_right_to_look = original_bp50_peak[0] - 5
    alternative_candidate_peaks_for35 = [peak for peak in highest_peaks_tup if peak[0] < furthest_right_to_look]
    best_guess = alternative_candidate_peaks_for35[-1]

    if best_guess[0] > original_bp35_peak[0]:
        log.write_to_log("looks like first peak of ladder might be hiding in the noise.")
        log.write_to_log("original_bp35_peak:" + str(original_bp35_peak[0]))
        log.write_to_log("alternative_candidate_peaks_for35:" + str(best_guess[0]))
        log.write_to_log("Will update first peak")
        return best_guess

    return original_bp35_peak



def remove_known_sus_ladder_peak(sixteen_peaks, highest_peaks_tup):
    sus_peaks = [x for x in sixteen_peaks if 2000 <= x[0] <= 2800]
    sus_peaks.sort(key=lambda x: x[1])
    num_peaks_to_remove = len(sus_peaks) - 3
    peaks_to_go = sus_peaks[0:num_peaks_to_remove]

    if num_peaks_to_remove > 0:
        print("peaks to go:" + str(peaks_to_go))
        print("sixteen_peaks:" + str(sixteen_peaks))
        # sixteen_peaks.remove(peaks_to_go)

        for i in range(0, num_peaks_to_remove):
            sixteen_peaks.remove(peaks_to_go[i])
            sixteen_peaks.append(highest_peaks_tup[i + 16])

    return sixteen_peaks


# best parameters for the ladder: (skinnier, cleaner spikes)
# kernel of 10, distance between peaks 2, min_peak_width 1;threshold_multiplier .35
# best parameters for the trace data: (fatter, messier spikes)
# kernel of 20, distance between peaks 20, min_peak_width 10;threshold_multiplier .5
def find_top_N_Peaks(signal_trace, num_peaks_to_find,
                                   kernel_size,
                                   min_distance_between_peaks,
                                   min_peak_width,
                                   threshold_multiplier, largest_first):  # .35
    # very basic smoothing
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_trace = np.convolve(signal_trace, kernel, mode='same')
    threshold = get_threshold_for_trace(smoothed_trace, threshold_multiplier)

    # https://plotly.com/python/peak-finding/
    indices = find_peaks(smoothed_trace, height=threshold, distance=min_distance_between_peaks, width=min_peak_width)[0]
    peak_heights_tup = [(x, smoothed_trace[x]) for x in indices]

    # sort, greatest peak height first. Keep the best N=num_peaks_to_find
    peak_heights_tup.sort(key=lambda x: x[1], reverse=True)
    highest_peaks_tup = peak_heights_tup[0:num_peaks_to_find]

    if (largest_first):
        return highest_peaks_tup, smoothed_trace, threshold
    else:
        highest_peaks_tup.sort(key=lambda x: x[0])
        return highest_peaks_tup, smoothed_trace, threshold

def get_threshold_for_trace(smoothed_trace, threshold_multiplier):
    # TODO, might be useful modification
    # calculate the threshold on range 2000:8000, skipping crazy peak
    # ladder_variance = np.var(smoothed_trace[2000:8000])
    # ladder_mode = mode(smoothed_trace[2000:8000])[0][0]
    ladder_variance = np.var(smoothed_trace)
    ladder_mode = mode(smoothed_trace)[0][0]
    ladder_sigma = np.sqrt(ladder_variance)
    threshold = ladder_mode + threshold_multiplier * ladder_sigma
    return threshold


def build_interpolation_based_on_ladder(run_folder, sixteen_peaks):
    # you are building a mapping from A -> B.
    # A = the peak positions in raw gel-travellijng space
    # B = the peak positions in bp length, as determined by the ladder.

    A = [x[0] for x in sixteen_peaks]  # get the peak position, not intensity
    B = [x[2] for x in sixteen_peaks]  # the LIZ 500 we baked in

    f1 = CubicSpline(A, B, bc_type='natural')

    # dont apply the spline where we dont have data!
    left_domain_limit = sixteen_peaks[0][0]
    right_domain_limit = sixteen_peaks[15][0]

    # test it looks good
    passed_monotonic_test, plotting_data_spline = per_file_visuals.plotMapping(run_folder, A, B, left_domain_limit,
                                                                               right_domain_limit,
                                                                               f1, "Spline")
    plotting_data_linear = False

    log.write_to_log("Ladder mapping is monotonic? " + str(passed_monotonic_test))

    if not passed_monotonic_test:
        log.write_warning_to_log("Ladder peaks are not well-spaced. Going to do linear interpolation instead of spine.")
        f1 = interp1d(A, B)
        log.write_to_log("Rebuilding and testing mapping... ")
        linear_mapping_plot_data = [A, B, left_domain_limit, right_domain_limit]
        passed_monotonic_test, plotting_data_linear = per_file_visuals.plotMapping(run_folder,
                                                                                   *linear_mapping_plot_data,
                                                                                   f1, "Linear")
        log.write_to_log("Ladder mapping is monotonic? " + str(passed_monotonic_test))

        if not passed_monotonic_test:
            return False

    return f1, left_domain_limit, right_domain_limit, plotting_data_spline, plotting_data_linear

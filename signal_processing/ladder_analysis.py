
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

    highest_peaks_tup, smoothed_trace, threshold = find_top_30_Peaks_largest_first(ladder_trace)

    #Option A
    #if we keep the largest:
    # of the remaining peaks, keep the greatest.
    highest_tup = highest_peaks_tup[0]

    # now, keep the best 15 starting from the right-most index.
    highest_peaks_tup.sort(key=lambda x: x[0], reverse=True)


    right_most = highest_peaks_tup[0:15]
    # put it together
    sixteen_peaks = [highest_tup] + right_most

    #Option B (option B currently seems to work better)
    #Highest 16 peaks, starting from the right
    sixteen_peaks = highest_peaks_tup[0:16]

    # want your ladder peaks leftmost on the left! not sorted by size
    sixteen_peaks.sort(key=lambda x: x[0])

    #Minor tweak, our vendor has a spurious peak that pops up between  2000 and 2800
    # SO - if we see 4 peaks between 2000 and 2800, and we should see ony 3,
    #throw out the smallest
    sixteen_peaks = remove_known_sus_ladder_peak(sixteen_peaks, highest_peaks_tup)
    sixteen_peaks.sort(key=lambda x: x[0])

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


def remove_known_sus_ladder_peak(sixteen_peaks, highest_peaks_tup):
    sus_peaks = [x for x in sixteen_peaks if 2000 <= x[0] <= 2800]
    sus_peaks.sort(key=lambda x: x[1])
    num_peaks_to_remove = len(sus_peaks) - 3
    peaks_to_go = sus_peaks[0:num_peaks_to_remove]

    if num_peaks_to_remove > 0:
        print("peaks to go:" + str(peaks_to_go))
        print("sixteen_peaks:" + str(sixteen_peaks))
        #sixteen_peaks.remove(peaks_to_go)

        for i in range(0, num_peaks_to_remove):
            sixteen_peaks.remove(peaks_to_go[i])
            sixteen_peaks.append(highest_peaks_tup[i+16])

    return sixteen_peaks


def find_top_30_Peaks_largest_first(ladder_trace):
    # very basic smoothing
    kernel_size = 20
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_trace = np.convolve(ladder_trace, kernel, mode='same')

    threshold = get_threshold_for_trace(smoothed_trace)

    # https://plotly.com/python/peak-finding/
    min_distance_between_peaks = 75
    min_peak_width = 10
    indices = find_peaks(smoothed_trace, height=threshold, distance=min_distance_between_peaks, width=min_peak_width)[0]
    peak_heights_tup = [(x, smoothed_trace[x]) for x in indices]
    # sort, greatest peak height first. Keep the best 30
    peak_heights_tup.sort(key=lambda x: x[1], reverse=True)
    highest_peaks_tup = peak_heights_tup[0:30]

    return highest_peaks_tup, smoothed_trace, threshold


def get_threshold_for_trace(smoothed_trace):
    # TODO, might be useful modification
    # calculate the threshold on range 2000:8000, skipping crazy peak
    # ladder_variance = np.var(smoothed_trace[2000:8000])
    # ladder_mode = mode(smoothed_trace[2000:8000])[0][0]
    ladder_variance = np.var(smoothed_trace)
    ladder_mode = mode(smoothed_trace)[0][0]
    ladder_sigma = np.sqrt(ladder_variance)
    threshold = ladder_mode + 0.35 * ladder_sigma
    return threshold


def build_interpolation_based_on_ladder(run_folder, sixteen_peaks):
    # you are building a mapping from A -> B.
    # A = the peak positions in raw gel-travellijng space
    # B = the peak positions in bp length, as determined by the ladder.

    A = [x[0] for x in sixteen_peaks]  # get the peak position, not intensity
    B = [x[2] for x in sixteen_peaks]  # the LIZ 500 we baked in

    f1 = CubicSpline(A, B, bc_type='natural')

    # dont appy the spline where we dont have data!
    # left_domain_limit = sixteen_peaks[0][0] - 100
    # right_domain_limit = sixteen_peaks[15][0] + 100
    left_domain_limit = sixteen_peaks[0][0]
    right_domain_limit = sixteen_peaks[15][0]

    #test it looks good
    passed_monotonic_test, plotting_data_spline = per_file_visuals.plotMapping(run_folder, A, B, left_domain_limit, right_domain_limit,
                                                                               f1, "Spline")
    plotting_data_linear=False

    log.write_to_log("Ladder mapping is monotonic? " + str(passed_monotonic_test))

    if not passed_monotonic_test:
        log.write_warning_to_log("Ladder peaks are not well-spaced. Going to do linear interpolation instead of spine.")
        f1 = interp1d(A, B)
        log.write_to_log("Rebuilding and testing mapping... ")
        linear_mapping_plot_data = [A, B, left_domain_limit, right_domain_limit]
        passed_monotonic_test, plotting_data_linear = per_file_visuals.plotMapping(run_folder, *linear_mapping_plot_data,
                                                                                   f1, "Linear")
        log.write_to_log("Ladder mapping is monotonic? " + str(passed_monotonic_test))

        if not passed_monotonic_test:
            return False

    return f1, left_domain_limit, right_domain_limit, plotting_data_spline, plotting_data_linear


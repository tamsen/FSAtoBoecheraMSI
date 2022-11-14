import numpy as np
from scipy.stats import mode
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from signal_processing import shared
from visualization import per_file_visuals
import log


GLOBAL_Liz500 = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500]

class Mapping_Info:

    def __init__(self, mapping_fxn, left_domain_limit, right_domain_limit,
                 mapping_plot_data_spline, mapping_plot_data_linear, results_are_questionable):

        self.mapping_fxn = mapping_fxn
        self.left_ladder_domain_limit = left_domain_limit
        self.right_ladder_domain_limit = right_domain_limit
        self.mapping_plot_data_spline = mapping_plot_data_spline
        self.mapping_plot_data_linear = mapping_plot_data_linear
        self.results_are_questionable = results_are_questionable

def getLadderPeaks(runFolder, runName, trace_data_dictionary):
    log.write_to_log("Reading through ladder trace for " + runName)
    # 'DATA105' is the ladder channel
    ladder_trace = trace_data_dictionary['DATA105']

    # parameters_for_right_side_of_ladder
    # parameters_for_right_side_of_ladder = shared.peak_calling_parameters(30, 20, 50, 10, .5)
    # parameters_for_left_side_of_ladder = shared.peak_calling_parameters(30, 20, 2, 1, .5)

    #made kernel smaller so ladder matches peak-smoothing alg
    parameters_for_right_side_of_ladder = shared.peak_calling_parameters(30, 10, 50, 10, .5, False)
    parameters_for_left_side_of_ladder = shared.peak_calling_parameters(30, 10, 2, 1, .5, False)

    highest_peaks_tup, smoothed_trace, threshold = find_top_N_Peaks(
        ladder_trace, parameters_for_right_side_of_ladder, True)

    highest_peaks_tup.sort(key=lambda x: x[0], reverse=True)

    if (len(highest_peaks_tup) < 18): #something went wrong here. Re-call the peaks

        new_theshold = highest_peaks_tup[0][1]*0.50
        parameters_for_right_side_of_ladder = shared.peak_calling_parameters(30, 10, 50, 10, .5, new_theshold)
        highest_peaks_tup, smoothed_trace, threshold = find_top_N_Peaks(
            ladder_trace, parameters_for_right_side_of_ladder, True)
        highest_peaks_tup.sort(key=lambda x: x[0], reverse=True)


    # now, keep the best 16 starting from the right-most index.
    highest_peaks_tup.sort(key=lambda x: x[0], reverse=True)
    sixteen_peaks = highest_peaks_tup[0:16]

    # want your ladder peaks leftmost on the left! not sorted by size
    sixteen_peaks.sort(key=lambda x: x[0])

    # Minor tweak, our vendor has a spurious peak that pops up between  2000 and 2800
    # SO - if we see 4 peaks between 2000 and 2800, and we should see ony 3,
    # throw out the smallest
    sixteen_peaks = remove_known_sus_ladder_peak_in_sus_range(sixteen_peaks, highest_peaks_tup, [2000,2800])
    sixteen_peaks.sort(key=lambda x: x[0])

    sixteen_peaks = remove_shorties_relative_to_siblings(sixteen_peaks, highest_peaks_tup)
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

    highest_peaks_tup, smoothed_trace, threshold = find_top_N_Peaks(
            ladder_trace, parameters_for_left_side_of_ladder, False)

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

def remove_shorties_relative_to_siblings(sixteen_peaks,highest_peaks_tup):

    expected_num_peaks=len(sixteen_peaks)
    indexes_to_check=range(4,len(sixteen_peaks)-3)
    peaks_to_go=[]
    significantly_shorter= 0.6

    for i in indexes_to_check:
        previous_peak_height=sixteen_peaks[i-1][1]
        ith_peak_height=sixteen_peaks[i][1]
        next_peak_height=sixteen_peaks[i+1][1]

        if (ith_peak_height < significantly_shorter*previous_peak_height) and \
            (ith_peak_height < significantly_shorter * next_peak_height):

            peaks_to_go.append(sixteen_peaks[i])

    num_peaks_to_remove = len(peaks_to_go)
    if num_peaks_to_remove > 0:

        for i in range(0, num_peaks_to_remove):

            sixteen_peaks.remove(peaks_to_go[i])
            sixteen_peaks.append(highest_peaks_tup[i + expected_num_peaks])

    return sixteen_peaks

def remove_known_sus_ladder_peak_in_sus_range(sixteen_peaks, highest_peaks_tup, sus_range):

    expected_num_peaks=len(sixteen_peaks)
    typical_sus_peak_BP100=sixteen_peaks[3]

    #sus_range=[2000,2800]
    #is sus peak in sus range?
    sus_peak_in_sus_range= sus_range[0] < typical_sus_peak_BP100[0] < sus_range[1]
    if not sus_peak_in_sus_range:
        return sixteen_peaks


    #is sus peak too close to its neightbors?
    # the distance between bp100 and bp139 is supposed to be significantly
    # bigger than
    # the distance between bp75 and bp100, and between bp139 and bp150.
    significantly_bigger= 1.6
    BP75=sixteen_peaks[2]
    BP139=sixteen_peaks[4]
    BP150=sixteen_peaks[5]
    dist75to100=typical_sus_peak_BP100[0]-BP75[0]
    dist100to139= BP139[0]- typical_sus_peak_BP100[0]
    dist139to50=BP150[0] - BP139[0]
    if (dist100to139 > significantly_bigger*dist75to100) and (dist100to139 > significantly_bigger*dist139to50):
        return sixteen_peaks

    sus_peaks = []
    sus_peak_index = []
    for i in range(0,expected_num_peaks):
        peak = sixteen_peaks[i]
        if sus_range[0] <= peak[0] <= sus_range[1] :
            sus_peaks.append(peak)
            sus_peak_index.append(i)

    sus_start=min(sus_peak_index)
    sus_end=max(sus_peak_index)
    peak_before_sus =sixteen_peaks[sus_start-1]
    peak_after_sus =sixteen_peaks[sus_end+1]
    height_of_shortest_nieghbor=min(peak_before_sus[1],peak_after_sus[1])

    sus_peaks.sort(key=lambda x: x[1])
    num_peaks_to_remove = len(sus_peaks) - 3
    peaks_to_go = sus_peaks[0:num_peaks_to_remove]

    if num_peaks_to_remove > 0:

        for i in range(0, num_peaks_to_remove):

            peak_to_go=peaks_to_go[i]
            if peak_to_go[1] < .95 * height_of_shortest_nieghbor:
                sixteen_peaks.remove(peaks_to_go[i])
                sixteen_peaks.append(highest_peaks_tup[i + expected_num_peaks])

    return sixteen_peaks

def find_top_N_Peaks(signal_trace, peak_calling_parameters,largest_first):

    # very basic smoothing
    kernel = np.ones(peak_calling_parameters.kernel_size) / peak_calling_parameters.kernel_size
    smoothed_trace = np.convolve(signal_trace, kernel, mode='same')

    if peak_calling_parameters.forced_threshold:
        threshold = peak_calling_parameters.forced_threshold
    elif peak_calling_parameters.threshold_multiplier:
        threshold = get_threshold_for_trace(smoothed_trace, peak_calling_parameters.threshold_multiplier)
    else:
        threshold = 0

    # https://plotly.com/python/peak-finding/
    indices = find_peaks(smoothed_trace, height=threshold,
                         distance=peak_calling_parameters.min_distance_between_peaks,
                         width=peak_calling_parameters.min_peak_width)[0]
    peak_heights_tup = [(x, smoothed_trace[x]) for x in indices]

    # sort, greatest peak height first. Keep the best N=num_peaks_to_find
    peak_heights_tup.sort(key=lambda x: x[1], reverse=True)
    highest_peaks_tup = peak_heights_tup[0:peak_calling_parameters.num_peaks_needed]

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

    # dont apply the spline where we dont have data!
    left_domain_limit = sixteen_peaks[0][0]
    right_domain_limit = sixteen_peaks[15][0]

    plotting_data_spline = False
    plotting_data_linear = False

    try:
        f1 = CubicSpline(A, B, bc_type='natural')

        # test it looks good
        passed_monotonic_test, plotting_data_spline = per_file_visuals.plotMapping(run_folder, A, B,
                                                                                   left_domain_limit,
                                                                                   right_domain_limit,
                                                                                   f1, "Spline")

        log.write_to_log("Ladder mapping is monotonic? " + str(passed_monotonic_test))

    except ValueError:
        log.write_to_log("Major spline fail: " + str(ValueError))
        passed_monotonic_test = False


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

    ladder_info = Mapping_Info(f1, left_domain_limit, right_domain_limit, plotting_data_spline, plotting_data_linear,
                               False)

    return ladder_info

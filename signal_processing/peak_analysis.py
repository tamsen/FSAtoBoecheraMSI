import statistics

import numpy as np


def peaks_to_raw_calls(peaks, trace_x_new, trace_y_new, threshold):
    peaks = [[round(peak[0], 1), peak[1]] for peak in peaks]
    return peaks


def peaks_to_filtered_calls(peaks, loci):

    #remove any distant peaks tht are 2/10 of the main peak
    iteratively_clean_short_peaks(peaks, 100, .2, 100, .2)

    step_width_left = 10
    step_proportion_left = 0.4
    step_width_right = 3
    step_proportion_right = 0.4

    #remove any close peaks that are 4/10 of the main peak
    iteratively_clean_short_peaks(peaks, step_width_left, step_proportion_left,
                                   step_width_right, step_proportion_right)


    typical_stutter = 3.5
    fine_stutter = 1.5

    # ----------- PS1 extra filtering --------------

    if loci == 'ICE3':
        peaks = stutter_fix(peaks, typical_stutter, take_run_maximum)

    if loci == 'BF20':  # most pain-in-the-ass loci
        merge_peaks_closer_than_this = fine_stutter  # g
        peaks = stutter_fix(peaks, merge_peaks_closer_than_this, take_run_maximum)

    # ----------- PS2 extra filtering --------------

    if loci in ['BF11', 'C8']:
        peaks = stutter_fix(peaks, typical_stutter, take_run_maximum)

    if loci in ['ICE14']:  # be careful here b/c we later remove peak at 210
        peaks = stutter_fix(peaks, fine_stutter, take_run_maximum)

    # ----------- PS3 extra filtering --------------

    if loci in ['BF9']:  # b
        peaks = stutter_fix(peaks, typical_stutter, take_run_maximum)

    if loci in ['E9']:  # b
        peaks = stutter_fix(peaks, typical_stutter, take_right_most)

    if loci == 'BF18':
        peaks = stutter_fix(peaks, typical_stutter, take_run_maximum)

    # ----------- PS4 extra filtering --------------

    if loci == 'BF3':  # second-most pain-in-the-ass loci
        merge_peaks_closer_than_this = 1
        peaks = stutter_fix(peaks, merge_peaks_closer_than_this, take_run_maximum)

    if loci == 'BF19':
        peaks = stutter_fix(peaks, typical_stutter, take_run_maximum)

    if loci == 'B6':
        peaks = stutter_fix(peaks, fine_stutter, take_run_maximum)

    # ----------- PS5 extra filtering --------------

    if loci == 'BF15':
        peaks = stutter_fix(peaks, typical_stutter, drop_right_most)

    if loci == 'Bdru266':
        peaks = stutter_fix(peaks, typical_stutter, take_run_maximum)

    return [[round(peak[0], 1), peak[1]] for peak in peaks]


def filter_by_range(peak_x, peak_y, expected_range):
    peaks_in_loci_range = []
    for i in range(0, len(peak_x)):
        x = peak_x[i]
        y = peak_y[i]
        if ((x >= expected_range[0]) and (x <= expected_range[1])):
            peaks_in_loci_range.append([x, y])

    peaks_in_loci_range.sort(key=lambda x: x[0])
    return peaks_in_loci_range


# Some times there is a teeny tiny peak to the left or rigth of a main peak.
# This code removes it.
def step_check(peaks, left_width, left_step_proportion,
               right_width, right_step_proportion):
    if len(peaks) == 0:
        return peaks

    xs = [p[0] for p in peaks]
    xs_diffs = np.diff(xs)
    ps_to_remove = []

    for i in range(0, len(peaks) - 1):

        d = xs_diffs[i]
        p_last = peaks[i]
        p_now = peaks[i + 1]
        if d <= left_width:
            if p_now[1] * left_step_proportion > p_last[1]:
                ps_to_remove.append(p_last)

    for i in range(0, len(peaks) - 1):

        d = xs_diffs[i]
        p_last = peaks[i]
        p_now = peaks[i + 1]
        if d <= right_width:
            if p_last[1] * right_step_proportion > p_now[1]:
                ps_to_remove.append(p_now)

    for p in ps_to_remove:

        if p in peaks:
            peaks.remove(p)

    return peaks


def iteratively_clean_short_peaks(peaks, step_width_left, step_proportion_left,
                                  step_width_right, step_proportion_right):

    i=0
    while True:

        oringinal_len = len(peaks)
        peaks = step_check(peaks, step_width_left, step_proportion_left,
                           step_width_right, step_proportion_right)

        new_len = len(peaks)
        #print("cleaned peaks: " + str(peaks))
        #print("iteration: " + str(i))
        i = i + 1

        if oringinal_len == new_len:
            break

    return peaks


def stutter_fix(peaks, how_close_is_too_close, f):
    if len(peaks) == 0:
        return peaks

    xs = [p[0] for p in peaks]
    xs_diffs = np.diff(xs)

    stutter_runs = [[]]
    peaks_in_stutter_runs = [[peaks[0]]]

    for i in range(0, len(peaks) - 1):

        d = xs_diffs[i]
        pi = peaks[i + 1]
        if d <= how_close_is_too_close:
            peaks_in_stutter_runs[-1].append(pi)
            stutter_runs[-1].append(d)

        if d > how_close_is_too_close:
            peaks_in_stutter_runs.append([pi])
            stutter_runs.append([d])

    best_peaks_in_a_run = []
    for run in peaks_in_stutter_runs:

        half_run_max = (take_run_maximum(run)[0])[1] * 0.5
        # if any peak is lest than 1/2 the max height in a run, kill it
        for p in run:
            if (p[1] < half_run_max):
                run.remove(p)

        best_peaks_in_a_run = best_peaks_in_a_run + f(run)

    return best_peaks_in_a_run


def take_run_maximum(run):
    if len(run) < 2:
        return run

    max_height_in_run = max([p[1] for p in run])
    for p in run:
        if p[1] == max_height_in_run:
            return [p]



def bf15_special(run):
    if len(run) < 2:
        return run

    min_height_in_run = min([p[1] for p in run])
    peaks_at_close_to_that_height = [p for p in run if min_height_in_run*(.95) < p[1] < min_height_in_run*(1.05)]

    p_replacement_x = statistics.mean( [p[0] for p in peaks_at_close_to_that_height])
    p_replacement_y = statistics.mean( [p[1] for p in peaks_at_close_to_that_height])

    for p in peaks_at_close_to_that_height:
        run.remove(p)

    if len(peaks_at_close_to_that_height) > 1:
        run.append([p_replacement_x,p_replacement_y])

    return run

def drop_lowest(run):
    if len(run) < 2:
        return run

    min_height_in_run = min([p[1] for p in run])
    to_drop = []
    for p in run:
        if p[1] == min_height_in_run:
            to_drop.append(p)

    for p in to_drop:
        run.remove(p)

    return run


def take_right_most(run):
    return [run[-1]]


def drop_right_most(run):
    if len(run) < 2:
        return (run)

    return run[0:-1]


def take_run_centroid(run):
    weight_sum = sum([p[1] for p in run])
    max_y = max([p[1] for p in run])
    x_weighted_avg = 0
    for p in run:
        x = p[0] * p[1]
        x_weighted_avg = x_weighted_avg + x
    return [[x_weighted_avg / weight_sum, max_y]]

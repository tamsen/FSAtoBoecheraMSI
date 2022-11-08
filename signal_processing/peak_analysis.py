import log
import numpy as np


def peaks_to_raw_calls(peaks, trace_x_new, trace_y_new, threshold):
    peaks = [[round(peak[0], 1), peak[1]] for peak in peaks]
    return  peaks


def peaks_to_filtered_calls(peaks, loci):
    left_step_width = 3
    left_step_proportion = 0.6
    peaks = left_step_check(peaks, left_step_width, left_step_proportion)

    # ----------- PS1 extra filtering --------------

    if loci == 'ICE3':
        merge_peaks_closer_than_this = 3.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this, take_run_centroid)

    if loci == 'BF20':
        # print("BF20")
        # print(str(peaks))
        merge_peaks_closer_than_this = 1.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this, take_run_maximum)

    # ----------- PS3 extra filtering --------------

    if loci == 'BF9':
        merge_peaks_closer_than_this = 4.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this, take_run_maximum)

    if loci == 'E9':
        merge_peaks_closer_than_this = 2.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this, take_right_most)

    # ----------- PS4 extra filtering --------------

    # none needed

    # ----------- PS5 extra filtering --------------

    if loci == 'BF15':
        merge_peaks_closer_than_this = 2
        peaks = stutter_check_2(peaks,  merge_peaks_closer_than_this, drop_right_most)

    if loci == 'Bdru266':
        merge_peaks_closer_than_this = 3
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this, drop_lowest)

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

def left_step_check(peaks, how_close_is_too_close, left_step_proportion):
    if len(peaks) == 0:
        return peaks

    xs = [p[0] for p in peaks]
    ys = [p[1] for p in peaks]
    print("xs: " + str(xs))
    print("ys: " + str(ys))
    xs_diffs = np.diff(xs)
    print("xs_diffs: " + str(xs_diffs))
    ps_to_remove = []

    for i in range(0, len(peaks) - 1):

        d = xs_diffs[i]
        p_last = peaks[i]
        p_now = peaks[i + 1]
        if d <= how_close_is_too_close:
            if p_now[1] * left_step_proportion > p_last[1]:
                ps_to_remove.append(p_last)

    for p in ps_to_remove:
        peaks.remove(p)

    return peaks


def stutter_check_2(peaks, how_close_is_too_close, f):
    if len(peaks) == 0:
        return peaks

    xs = [p[0] for p in peaks]
    ys = [p[1] for p in peaks]
    print("xs: " + str(xs))
    print("ys: " + str(ys))
    xs_diffs = np.diff(xs)
    print("xs_diffs: " + str(xs_diffs))

    stutter_runs = [[]]
    peaks_in_stutter_runs = [[peaks[0]]]

    for i in range(0, len(peaks) - 1):

        d = xs_diffs[i]
        pi = peaks[i + 1]
        print("d:" + str(d))
        if d <= how_close_is_too_close:
            peaks_in_stutter_runs[-1].append(pi)
            stutter_runs[-1].append(d)

        if d > how_close_is_too_close:
            peaks_in_stutter_runs.append([pi])
            stutter_runs.append([d])

    best_peaks_in_a_run = []
    print("stutter_runs: " + str(stutter_runs))
    print("peaks_in_stutter_runs: ")
    for run in peaks_in_stutter_runs:
        print("run:  " + str(run))

        half_run_max = (take_run_maximum(run)[0])[1] * 0.5
        # if any peak is lest than 1/2 the max height in a run, kill it
        for p in run:
            if (p[1] < half_run_max):
                run.remove(p)

        best_peaks_in_a_run = best_peaks_in_a_run + f(run)

    print("best peaks: " + str(best_peaks_in_a_run))
    return best_peaks_in_a_run


def take_run_maximum(run):
    print(str(run))
    max_height_in_run = max([p[1] for p in run])
    for p in run:
        if p[1] == max_height_in_run:
            return [p]


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
    return  [ run[-1] ]

def drop_right_most(run):

    if len(run) < 2 :
        return(run)

    return run[0:-1]


def take_run_centroid(run):
    weight_sum = sum([p[1] for p in run])
    max_y = max([p[1] for p in run])
    x_weighted_avg = 0
    for p in run:
        x = p[0] * p[1]
        x_weighted_avg = x_weighted_avg + x
    return [ [x_weighted_avg / weight_sum, max_y] ]

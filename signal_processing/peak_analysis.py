import statistics

import numpy as np


def peaks_to_raw_calls(peaks):
    peaks = [[round(peak[0], 1), peak[1], peak[2]] for peak in peaks]
    return peaks


def peaks_to_filtered_calls(peaks, loci, trace_x_new, trace_y_new, threshold):


    #if loci=="BF3":
    #    print("break here")

    #brutally remove any distant peaks, and get gentler, closer to the action
    brutally_clean_short_peaks(peaks, .1)
    iteratively_clean_short_peaks(peaks, 100, .2, 100, .2)
    iteratively_clean_short_peaks(peaks, 10, 0.4, 3,0.4)

    typical_stutter = 3.5
    fine_stutter = 1.5

    # ----------- PS1 extra filtering --------------

    if loci == 'ICE3':
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, take_run_maximum)

    if loci == 'BF20':  # most pain-in-the-ass loci
        merge_peaks_closer_than_this = fine_stutter  # g
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, merge_peaks_closer_than_this, take_run_maximum)

    # ----------- PS2 extra filtering --------------

    if loci in ['BF11', 'C8']:
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, take_run_maximum)

    if loci in ['ICE14']:  # be careful here b/c we later remove peak at 210
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, fine_stutter, take_run_maximum)

    # ----------- PS3 extra filtering --------------

    if loci in ['BF9']:  # b
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold,
                            typical_stutter, bf9_special)
        num_peaks=len(peaks)
        if num_peaks > 3:
            peaks=peaks[num_peaks-4:num_peaks]

    if loci in ['E9']:  # b
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, take_right_most)

    if loci == 'BF18':
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, take_run_maximum)

    # ----------- PS4 extra filtering --------------

    if loci == 'BF3':  # second-most pain-in-the-ass loci
        merge_peaks_closer_than_this = 1
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold,
                            merge_peaks_closer_than_this, take_run_maximum)

    if loci == 'BF19':
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, take_run_maximum)

    if loci == 'B6':
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, fine_stutter, take_run_maximum)

    # ----------- PS5 extra filtering --------------

    if loci == 'BF15':
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, drop_right_most)

    if loci == 'Bdru266':
        peaks = stutter_fix(peaks, trace_x_new, trace_y_new, threshold, typical_stutter, take_run_maximum)

    return [[round(peak[0], 1), peak[1]] for peak in peaks]


def get_min_in_trace_section(x1,x2,trace_x_new, trace_y_new):

    ys=[]
    for i in range(0,len(trace_x_new)):
        x=trace_x_new[i]
        y=trace_y_new[i]
        if x1 < x < x2:
            ys.append(y)
    return min(ys)

def bf9_special(run,trace_x_new, trace_y_new, threshold):

    if len(run) < 2:
        return run

    peaks_to_remove = []
    for i in range (0,len(run)-1):

         this_peak=run[i]
         next_peak=run[i+1]
         this_peak_x=this_peak[0]
         next_peak_x=next_peak[0]

         min_y_between_peaks = get_min_in_trace_section(this_peak_x,next_peak_x,trace_x_new, trace_y_new)
         dip_between_peaks= min_y_between_peaks < threshold

         if not dip_between_peaks:

             if this_peak[1] <= next_peak[1]:
                 peaks_to_remove.append(this_peak)
             else:
                 peaks_to_remove.append(next_peak)

    for peak in peaks_to_remove:
        if peak in run:
            run.remove(peak)

    #keep the right-most two
    num_left=len(run)
    return run[num_left-2:num_left]

def filter_by_range(peak_x, peak_y, dip_between_peaks, expected_range):
    peaks_in_loci_range = []
    for i in range(0, len(peak_x)):
        x = peak_x[i]
        y = peak_y[i]
        dip = dip_between_peaks[i]
        if ((x >= expected_range[0]) and (x <= expected_range[1])):
            peaks_in_loci_range.append([x, y, dip])

    peaks_in_loci_range.sort(key=lambda x: x[0])
    return peaks_in_loci_range


# Some times there is a teeny tiny peak to the left or right of a main peak.
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

def brutally_clean_short_peaks(peaks, step_proportion):

    if len(peaks) < 2:
        return peaks

    max_peak = max([p[1] for p in peaks])
    peaks_to_go=  [p for p in peaks if p[1] < step_proportion*max_peak]

    for p in peaks_to_go:
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


def stutter_fix(peaks,trace_x_new, trace_y_new, threshold, how_close_is_too_close, f):
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

        #tjd - should I put this back? or, ony remove if there is a significant dip?
        #half_run_max = (take_run_maximum(run,trace_x_new, trace_y_new, threshold)[0])[1] * 0.5
        # if any peak is lest than 1/2 the max height in a run, kill it
        #for p in run:
        #    if (p[1] < half_run_max):
        #        run.remove(p)

        best_peaks_in_a_run = best_peaks_in_a_run + f(run, trace_x_new, trace_y_new, threshold)

    return best_peaks_in_a_run


def take_run_maximum(run, trace_x_new, trace_y_new, threshold):
    if len(run) < 2:
        return run

    max_height_in_run = max([p[1] for p in run])
    for p in run:
        if p[1] == max_height_in_run:
            return [p]



def drop_lowest(run,trace_x_new, trace_y_new, threshold):
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


def take_right_most(run, trace_x_new, trace_y_new, threshold):
    return [run[-1]]


def drop_right_most(run, trace_x_new, trace_y_new, threshold):
    if len(run) < 2:
        return (run)

    return run[0:-1]


def take_run_centroid(run, trace_x_new, trace_y_new, threshold):
    weight_sum = sum([p[1] for p in run])
    max_y = max([p[1] for p in run])
    x_weighted_avg = 0
    for p in run:
        x = p[0] * p[1]
        x_weighted_avg = x_weighted_avg + x
    return [[x_weighted_avg / weight_sum, max_y]]

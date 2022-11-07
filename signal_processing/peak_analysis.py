import log
import numpy as np

def peaks_to_raw_calls(peaks, trace_x_new, trace_y_new, threshold):

    MSI_calls = [[round(peak[0],1), peak[1]] for peak in peaks]

    return MSI_calls

def peaks_to_filtered_calls(peaks, loci, trace_x_new, trace_y_new, threshold):

    left_step_width=3
    left_step_proportion=0.6
    peaks = left_step_check(peaks, left_step_width,left_step_proportion)

    #----------- PS1 extra filtering --------------

    if loci == 'ICE3':
        merge_peaks_closer_than_this = 3.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this,take_run_centroid)

    if loci == 'BF20':
        #print("BF20")
        #print(str(peaks))
        merge_peaks_closer_than_this = 1.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this,take_run_maximum)

    #----------- PS3 extra filtering --------------

    if loci == 'BF9':
        merge_peaks_closer_than_this = 4.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this,take_run_maximum)

    if loci == 'E9':
        merge_peaks_closer_than_this = 2.5
        peaks = stutter_check_2(peaks, merge_peaks_closer_than_this,take_right_most)


    return [[round(peak[0],1), peak[1]] for peak in peaks]

def filter_by_range(peak_x, peak_y, expected_range):

    peaks_in_loci_range=[]
    for i in range(0, len(peak_x)):
            x = peak_x[i]
            y = peak_y[i]
            if  ((x >= expected_range[0]) and (x <= expected_range[1])):

                peaks_in_loci_range.append([x,y])

    peaks_in_loci_range.sort(key=lambda x: x[0])
    return peaks_in_loci_range


def find_minimum_between_peaks(Peak1, Peak2, trace_x_new, trace_y_new):

    lowest_y_value=max(trace_y_new)
    x_for_lowest_y_value=0

    for i in range(0,len(trace_x_new)):
        x= trace_x_new[i]
        if (x > Peak1[0]) and (x < Peak2[0]):
            y = trace_y_new[i]
            if (y <= lowest_y_value):
                lowest_y_value =y
                x_for_lowest_y_value = x

    log.write_to_log("PeakA:" + str(Peak1))
    log.write_to_log("PeakB:" + str(Peak2))
    log.write_to_log("Lowest Y:" + str(lowest_y_value))

    return [ x_for_lowest_y_value,lowest_y_value]


def there_is_never_a_drop_between_peaks(Peak1, Peak2, trace_x_new, trace_y_new, threshold,
                                        required_drop_fraction):

    [x_for_lowest_y_value, lowest_y_value] = find_minimum_between_peaks(Peak1, Peak2, trace_x_new, trace_y_new)

    drop_distance= abs( min( Peak1[1],Peak2[1]) - lowest_y_value )
    distance_between_peak_and_threshold =  abs( min( Peak1[1],Peak2[1]) - threshold )

    log.write_to_log("checking for drop between " + str(Peak1[0]) + " and " + str(Peak2[0]))

    log.write_to_log("drop_distance " + str(drop_distance))
    log.write_to_log("required_drop_fraction " + str(required_drop_fraction))
    log.write_to_log("distance_between_peak_and_threshold " + str(distance_between_peak_and_threshold))
    log.write_to_log("required_drop_fraction * distance_between_peak_and_threshold " + str(
        required_drop_fraction * distance_between_peak_and_threshold))
    if drop_distance < required_drop_fraction* distance_between_peak_and_threshold:
        log.write_to_log("No drop between peaks observed. They are probably the same peak.")
        return True
    else:
        log.write_to_log("Drop observed between peaks. They are probably true different peaks.")
        return False

def find_threshold_crossings(trace_x_new, trace_y_new, threshold):

    above_threshold_intervals=[]
    above_threshold=False
    interval_start=0;

    for i in range(0,trace_y_new):

        #if we start above the threshold, and pass below, close off the interval
        if above_threshold and (trace_y_new[i] <= threshold):
            above_threshold_intervals.append([interval_start,trace_x_new[i]])
            above_threshold = False

        #if we start below the threshold, and pass above, open the interval
        if (not above_threshold) and (trace_y_new[i+1] < threshold):
            interval_start = trace_x_new[i]
            above_threshold = True

    return above_threshold_intervals

def left_step_check(peaks, how_close_is_too_close,left_step_proportion):

    if len(peaks) == 0:
        return peaks

    xs=[p[0] for p in peaks]
    ys=[p[1] for p in peaks]
    print("xs: " + str(xs))
    print("ys: " + str(ys))
    xs_diffs=np.diff(xs)
    print ("xs_diffs: " + str(xs_diffs))
    ps_to_remove=[]

    for i in range(0,len(peaks)-1):

        d=xs_diffs[i]
        p_last= peaks[i]
        p_now= peaks[i+1]
        if d <= how_close_is_too_close:
            if p_now[1]*left_step_proportion > p_last[1]:
                ps_to_remove.append(p_last)

    for p in ps_to_remove:
        peaks.remove(p)

    return peaks

def stutter_check_2(peaks, how_close_is_too_close, f):

    if len(peaks) == 0:
        return peaks

    xs=[p[0] for p in peaks]
    ys=[p[1] for p in peaks]
    print("xs: " + str(xs))
    print("ys: " + str(ys))
    xs_diffs=np.diff(xs)
    print ("xs_diffs: " + str(xs_diffs))

    stutter_runs=[[]]
    peaks_in_stutter_runs=[[peaks[0]]]

    for i in range(0,len(peaks)-1):

        d=xs_diffs[i]
        pi= peaks[i+1]
        print("d:" + str(d))
        if d <= how_close_is_too_close:
            peaks_in_stutter_runs[-1].append(pi)
            stutter_runs[-1].append(d)

        if d > how_close_is_too_close:
            peaks_in_stutter_runs.append([pi])
            stutter_runs.append([d])

    best_peaks_in_a_run=[]
    print("stutter_runs: " + str(stutter_runs))
    print("peaks_in_stutter_runs: " )
    for run in peaks_in_stutter_runs:
        print("run:  " + str(run))

        half_run_max=take_run_maximum(run)[1]*0.5
        #if any peak is lest than 1/2 the max height in a run, kill it
        for p in run:
            if(p[1]<half_run_max):
                run.remove(p)

        best_peaks_in_a_run.append(f(run))

    print("stutter_run_maximun: " + str(best_peaks_in_a_run))
    return best_peaks_in_a_run

def take_run_maximum(run):
    max_height_in_run = max([p[1] for p in run])
    for p in run:
        if p[1] == max_height_in_run:
            return p

def take_right_most(run):
    return run[-1]

def take_run_centroid(run):

    weight_sum=sum([p[1] for p in run])
    max_y=max([p[1] for p in run])
    x_weighted_avg=0
    for p in run:
        x= p[0]  * p[1]
        x_weighted_avg = x_weighted_avg+x
    return [x_weighted_avg / weight_sum, max_y]


def insist_on_drops_between_peaks(peaks, trace_x_new, trace_y_new, threshold, required_drop_fraction):

    if (len(peaks)) <= 1:
        consolidated_peaks=peaks
    else:
        consolidated_peaks = []
    i = 0
    while i < (len(peaks)-1):

        xi0=peaks[i][0]
        yi0=peaks[i][1]
        xi1=peaks[i+1][0]
        yi1=peaks[i+1][1]
        log.write_to_log("checking peak " + str(xi0) + "has a dip between itself and other peaks")
        if there_is_never_a_drop_between_peaks(peaks[i], peaks[i + 1], trace_x_new, trace_y_new, threshold,
                                               required_drop_fraction):
            if yi1 > yi0:

                if (i+1 == len(peaks) -1): #we are on the last one, and its the biggest
                    consolidated_peaks.append( [ xi1, yi1])
                i = i + 1
            else:
                consolidated_peaks.append( [ xi0, yi0])
                i = i + 1

        else:
                consolidated_peaks.append(peaks[i])
                i = i +1
                if i == len(peaks)-1:
                    consolidated_peaks.append(peaks[i])

    return consolidated_peaks


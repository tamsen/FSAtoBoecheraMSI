import log

def peaks_to_msi_calls(peaks, trace_x_new, trace_y_new, threshold):

    #parameters
    merge_peaks_closer_than_this  = 1.5
    required_drop_fraction = 0.5

    #peaks = consolidate(peaks, merge_peaks_closer_than_this)
    #peaks = insist_on_drops_between_peaks(peaks, trace_x_new, trace_y_new, threshold, required_drop_fraction)
    #skipping stutter check, dont seem to need it?
    MSI_calls = [[round(peak[0],1), peak[1]] for peak in peaks]

    return MSI_calls

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

def consolidate(peaks, how_close_is_too_close):

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
        log.write_to_log("checking peaks " + str(peaks[i]) + "and" + str(peaks[i + 1]) + "for possible consolidation")
        if (xi1 - xi0) < how_close_is_too_close:
            consolidated_peaks.append( [ (xi1 + xi0) / 2.0 , (yi1 + yi0) / 2.0] )
            log.write_to_log("Consolidating " + str(peaks[i]) + "and" + str(peaks[i + 1]))
            i = i +2
        else:
            log.write_to_log("No consolidation required.")
            consolidated_peaks.append(peaks[i])
            i = i +1
            if i == len(peaks)-1:
                consolidated_peaks.append(peaks[i])

    return consolidated_peaks


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



def check_for_stutter(peaks):

    stutter_size=2
    stutter_peak_indexes=[]
    for i in range(0,len(peaks)-2):

        three_xs= [peak[0] for peak in peaks[i:i+3]]
        three_ys= [peak[1] for peak in peaks[i:i+3]]

        if (three_xs[1]-three_xs[0] <= stutter_size) \
           and (three_xs[2]-three_xs[1]  <= stutter_size):

            print("Peaks frequently 2bp apart detected..")
            print("peak_xs[i:i+3]:" + str(three_xs))
            print("peak_ys[i:i+3]:" + str(three_ys))

            #if they neatly have decreasing magnitude, we accept them.
            # Otherwise, possible stutter run

            if (three_ys[1] - three_ys[0] >= 0) \
                    or (three_ys[2] - three_ys[1] >= 0):

                #if we got hear, either y1 or y2 is the highest.
                if (three_ys[2] - three_ys[1] >= 0): #y2 is highest, so y0 and y1 are not real
                    stutter_peak_indexes = stutter_peak_indexes + [i,i+1]
                else: #y1 is highest, so y0 and y2 are not real
                    stutter_peak_indexes = stutter_peak_indexes + [i+1,i+2]

                log.write_to_log("Amplitude patterns doe not look like real alleles.")
                log.write_to_log("peaks[i:i+3]:" + str(three_ys))

    #keep unique stutter indexes
    stutter_peak_indexes = set(stutter_peak_indexes)
    de_stuttered_peaks=[]
    for i in range(0,len(peaks)):
        if i not in  stutter_peak_indexes:
            de_stuttered_peaks.append(peaks[i])

    return de_stuttered_peaks
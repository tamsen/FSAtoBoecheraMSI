
def PeaksToMsiCalls(peaks, trace_x_new, trace_y_new):

    peaks = Consolidate(peaks)
    peaks = CheckForStutter(peaks)
    MSI_calls = [[round(peak[0]), peak[1]] for peak in peaks]

    return MSI_calls

def FilterByRange(peak_x, peak_y, expected_range):

    peaks_in_loci_range=[]
    for i in range(0, len(peak_x)):
            x = peak_x[i]
            y = peak_y[i]
            if  ((x >= expected_range[0]) and (x <= expected_range[1])):

                peaks_in_loci_range.append([x,y])

    peaks_in_loci_range.sort(key=lambda x: x[0])
    return peaks_in_loci_range

def FindThresholdCrossings(peaks,trace_x_new, trace_y_new,threshold):

    crossings_up=[]
    crossings_down=[]
    above_threshold_intervals=[]
    for i in range(0,trace_y_new-1):

        if (trace_y_new[i] <= threshold) and (trace_y_new[i+1] > threshold):
            crossings_up.append(trace_x_new[i])
        if (trace_y_new[i] >= threshold) and (trace_y_new[i+1] < threshold):
            crossings_down.append(trace_x_new[i])

def Consolidate(peaks):

    how_close_is_too_close = 1.2
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

        if (xi1 - xi0) < how_close_is_too_close:
            consolidated_peaks.append( [ (xi1 + xi0) / 2.0 , (yi1 + yi0) / 2.0] )
            i = i +2
        else:
                consolidated_peaks.append(peaks[i])
                i = i +1
                if i == len(peaks)-1:
                    consolidated_peaks.append(peaks[i])

    #MSI_calls = [[round(peak[0]),peak[1]] for peak in consolidated_peaks]
    MSI_calls = [[peak[0], peak[1]] for peak in consolidated_peaks]

    return MSI_calls


def CheckForStutter(peaks):

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

                print("Amplitude patterns doe not look like real alleles.")
                print("peaks[i:i+3]:" + str(three_ys))

    #keep unique stutter indexes
    stutter_peak_indexes = set(stutter_peak_indexes)
    de_stuttered_peaks=[]
    for i in range(0,len(peaks)):
        if i not in  stutter_peak_indexes:
            de_stuttered_peaks.append(peaks[i])

    return de_stuttered_peaks
import log
from signal_processing import peak_analysis


def make_adjustments(peaks, loci):
    if loci == 'B6':
        return peaks

    merge_peaks_closer_than_this = 2.5
    required_drop_fraction = 0.5

    #fw1379 has issue with PS5, there is a bug in stutter check
    #peaks = peak_analysis.check_for_stutter(peaks)
    #but proably nees some code so for some peaks, if they are within 2bp, take the leading (smaller distance)
    # #or try the larger peak.

    peaks = consolidate(peaks, merge_peaks_closer_than_this)

    return peaks


def consolidate(peaks, how_close_is_too_close):
    if (len(peaks)) <= 1:
        consolidated_peaks = peaks
    else:
        consolidated_peaks = []
    i = 0
    while i < (len(peaks) - 1):

        xi0 = peaks[i][0]
        yi0 = peaks[i][1]
        xi1 = peaks[i + 1][0]
        yi1 = peaks[i + 1][1]
        log.write_to_log("checking peaks " + str(peaks[i]) + "and" + str(peaks[i + 1]) + "for possible consolidation")
        if (xi1 - xi0) < how_close_is_too_close:

            if yi1 > yi0:
                consolidated_peaks.append([xi1, yi1])
            else:
                consolidated_peaks.append([xi0, yi0])
            #consolidated_peaks.append([(xi1 + xi0) / 2.0, (yi1 + yi0) / 2.0])
            log.write_to_log("Consolidating " + str(peaks[i]) + "and" + str(peaks[i + 1]))
            i = i + 2
        else:
            log.write_to_log("No consolidation required.")
            consolidated_peaks.append(peaks[i])
            i = i + 1
            if i == len(peaks) - 1:
                consolidated_peaks.append(peaks[i])

    return consolidated_peaks

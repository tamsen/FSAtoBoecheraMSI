import log
from signal_processing import peak_analysis


def make_adjustments(peaks, loci):

    if loci == 'B6':
        return [[round(p[0], 0), p[1]] for p in peaks]

    #-------------PS1--------------

    if loci == 'ICE3':

        adjusted_peaks=[]
        for p in peaks:

            p0=p[0]
            if p0 < 100:
                p0=p[0]-3.0
            else:
                p0=p[0]-2.0

            adjusted_peaks.append([p0,p[1]])

        peaks = adjusted_peaks

    if loci == 'BF20':
        [peaks.remove(p) for p in peaks if 203 < p[0] < 204.5]
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    #-------------PS2--------------

    if loci == 'ICE14':
        [peaks.remove(p) for p in peaks if 209.5 < p[0] < 210.5]
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    if loci == 'BF11':
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]

    #-------------PS3--------------

    if loci == 'BF9':
        [peaks.remove(p) for p in peaks if 120 < p[0] < 122]

    if loci == 'E9':
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    return [[int(round(p[0], 0)), p[1]] for p in peaks]

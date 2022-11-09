import log
from signal_processing import peak_analysis

#removes "spurious" peaks and aligns calls to MW original binning
def make_adjustments(peaks, loci):

    #-------------PS1--------------

    if loci == 'ICE3':

        adjusted_peaks=[]
        for p in peaks:

            p0=p[0]
            if p0 < 100:
                p0=p[0]-1.7#3.0
            else:
                p0=p[0]-1.7 #2.0

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
        [peaks.remove(p) for p in peaks if 95 < p[0] < 96]
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]

    #-------------PS3--------------

    if loci == 'BF9':
        [peaks.remove(p) for p in peaks if 120 < p[0] < 122]
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]

    if loci == 'E9':
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]


    #-------------PS4--------------

    if loci in ['BF3','BF19']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]

    if loci in ['B6']:
        [peaks.remove(p) for p in peaks if 308 < p[0] < 310]
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    #-------------PS5--------------

    if loci in ['Bdru266']:
        [peaks.remove(p) for p in peaks if 80 < p[0] < 100]

    if loci in ['A3']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]


    return [[int(round(p[0], 0)), p[1]] for p in peaks]

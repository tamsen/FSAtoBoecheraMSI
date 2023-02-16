
#removes "spurious" peaks and aligns calls to MW original binning
def make_adjustments(peaks, loci):

    #-------------PS1--------------

    if loci == 'ICE3': #b
        peaks = [[p[0] -1.7, p[1]] for p in peaks]

    if loci == 'BF20': #g
        peaks = [[p[0] - 0.3, p[1]] for p in peaks]

    if loci == 'A1': #b
        peaks = [[p[0] -0.5, p[1]] for p in peaks]

    #-------------PS2--------------

    if loci == 'ICE14': #g
        peaks = [[p[0] - 0.0, p[1]] for p in peaks]

    if loci == 'BF11': #b
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]

    if loci == 'C8':
        peaks = [[p[0] - 0.6, p[1]] for p in peaks]

    #-------------PS3--------------

    if loci == 'BF9': #b
        peaks = [[p[0] - 6.0, p[1]] for p in peaks]

    if loci == 'BF18': #g
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]


    if loci == 'BF18': #g
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]

    #-------------PS4--------------

    if loci in ['BF3']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]

    if loci in ['BF19']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]

    if loci in ['B6']:
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    #-------------PS5--------------

    if loci in ['Bdru266']:
        [peaks.remove(p) for p in peaks if p[0] < 100]

    if loci in ['A3']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]

    if loci in ['BF15']:
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    return [[int(round(p[0], 0)), p[1]] for p in peaks]


#removes "spurious" peaks and aligns calls to MW original binning
def make_adjustments(peaks, loci):

    #-------------PS1--------------

    if loci == 'ICE3': #b
        peaks = [[p[0] -0.5, p[1]] for p in peaks]

        for p in peaks:
            if p[0] > 131:
               p[0] = p[0]-1.5

    if loci == 'BF20': #b

        peaks = [[p[0] +0.1, p[1]] for p in peaks]


    #-------------PS2--------------

    if loci == 'BF11': #b
        [peaks.remove(p) for p in peaks if 94 < p[0] < 98]
        peaks = [[p[0] + 1.5, p[1]] for p in peaks]

    if loci == 'ICE14': #g
        [peaks.remove(p) for p in peaks if 209.5 < p[0] < 210.5] #MW determined this to be false peak

    if loci == 'C8':
        [peaks.remove(p) for p in peaks if 94 < p[0] < 98]
        peaks = [[p[0] - 0.0, p[1]] for p in peaks]

    #-------------PS3--------------

    if loci == 'BF9': #b

        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    if loci == 'BF18': #g
        peaks = [[p[0] + 1.5, p[1]] for p in peaks]



    #-------------PS4--------------

    if loci in ['BF3']:
        peaks = [[p[0] - 0.0, p[1]] for p in peaks]

    if loci in ['BF19']:
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]
        [peaks.remove(p) for p in peaks if p[0] < 120]  # known false peaks, noise

    if loci in ['B6']:
        [peaks.remove(p) for p in peaks if 308 < p[0] < 310]
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    #-------------PS5--------------

    if loci in ['Bdru266']:
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    if loci in ['A3']:
        peaks = [[p[0] - 0.0, p[1]] for p in peaks]

    if loci in ['BF15']:
        peaks = [[p[0] + 0.0, p[1]] for p in peaks]

    return [[int(round(p[0], 0)), p[1]] for p in peaks]

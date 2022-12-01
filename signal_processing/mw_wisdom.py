
#removes "spurious" peaks and aligns calls to MW original binning
def make_adjustments(peaks, loci):

    #-------------PS1--------------

    if loci == 'ICE3': #b

        peaks = [[p[0] -1.7, p[1]] for p in peaks]
        # In original meeting, discussing GeneMarker results, we saw p > 100 had to be shifted -3
        # but using my script seems to require less shift. (meeting 2022 Aug 09 )
        # MW proposed capping the allele at 131, when trying to match the database.
        # (Because this allele is constantly marching up in known lemmonii populations) but again,
        # I dont think this cap is needed. TBD.
        #adjusted_peaks=[]
        #for p in peaks:

        #    p0=p[0]
        #    if p0 < 131:
        #        p0 = p[0]-1.7
        #    else:
        #        p0 = 131

        #    adjusted_peaks.append([p0,p[1]])

        #peaks = adjusted_peaks

    if loci == 'BF20': #g
        [peaks.remove(p) for p in peaks if 203 < p[0] < 206]  # MW determined this to be false peak
        [peaks.remove(p) for p in peaks if 235 < p[0] ]  #cross contamination

        peaks = [[p[0] - 0.0, p[1]] for p in peaks]

    if loci == 'A1': #b
        peaks = [[p[0] -0.5, p[1]] for p in peaks]

    #-------------PS2--------------

    if loci == 'ICE14': #g
        [peaks.remove(p) for p in peaks if 209.5 < p[0] < 210.5] #MW determined this to be false peak

    if loci == 'BF11': #b
        [peaks.remove(p) for p in peaks if 95 < p[0] < 96]
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]

    #-------------PS3--------------

    if loci == 'BF9': #b
        [peaks.remove(p) for p in peaks if 120 < p[0] < 122] # known false peaks, noise
        [peaks.remove(p) for p in peaks if 140 < p[0] < 150] # known false peaks, noise
        #[peaks.remove(p) for p in peaks if 86 < p[0] < 88] # known false peaks, noise
        peaks = [[p[0] - 6.0, p[1]] for p in peaks]

        if (len(peaks) > 3):
            [peaks.remove(p) for p in peaks if 86 < p[0] < 89]  # known false peaks, noise
            [peaks.remove(p) for p in peaks if 72 < p[0] < 78]  # known false peaks, noise

        if (len(peaks) > 3):
            peaks.remove(peaks[0]) #if all else fails, usually the last peak is the false one

    if loci == 'BF18': #g
        #peaks = [[p[0] - 1, p[1]] for p in peaks]
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]
        #peaks = [[p[0] - 0.0, p[1]] for p in peaks]

    if loci == 'BF18': #g
        #peaks = [[p[0] - 1, p[1]] for p in peaks]
        peaks = [[p[0] - 0.5, p[1]] for p in peaks]
        #peaks = [[p[0] - 0.0, p[1]] for p in peaks]
    #-------------PS4--------------

    if loci in ['BF3']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]

    if loci in ['BF19']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]
        [peaks.remove(p) for p in peaks if p[0] < 120]  # known false peaks, noise

    if loci in ['B6']:
        [peaks.remove(p) for p in peaks if 308 < p[0] < 310]
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    #-------------PS5--------------

    if loci in ['Bdru266']:
        [peaks.remove(p) for p in peaks if p[0] < 100]

    if loci in ['A3']:
        peaks = [[p[0] - 1.0, p[1]] for p in peaks]

    if loci in ['BF15']:
        peaks = [[p[0] + 0.5, p[1]] for p in peaks]

    return [[int(round(p[0], 0)), p[1]] for p in peaks]

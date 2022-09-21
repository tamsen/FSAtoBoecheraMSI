

def FilterByRange(peak_x, peak_y, expected_range):

    peaks_in_loci_range=[]
    for i in range(0, len(peak_x)):
            x = peak_x[i]
            y = peak_y[i]
            if  ((x >= expected_range[0]) and (x <= expected_range[1])):

                peaks_in_loci_range.append([x,y])

    peaks_in_loci_range.sort(key=lambda x: x[0])

def Consolidate(peaks):

    consolidated_peaks=[]
    i = 0
    while i < (len(peaks)-1):

        xi0=peaks[0][i]
        yi0=peaks[1][i]
        xi1=peaks[0][i+1]
        yi1=peaks[1][i+1]

        if (xi1 - xi0) < 1:
            consolidated_peaks.append( [ xi1 - xi0 / 2.0 , yi1 - yi0 / 2.0] )
            i = i +2
        else:
                consolidated_peaks.append(peaks[i])
                i = i +1
                if i == len(peaks)-1:
                    consolidated_peaks.append(peaks[i])

    #Check for stutter - do we have a sequence of peaks 2bp apart?
    #Is the height always increasing over the stutter?
    #if yes, skip the lower (left) peaks, and take the highest (right-most)



    MSI_calls = [round(peak) for peak in consolidated_peaks]
    MSI_calls = [ round(peak,1) for peak in consolidated_peaks]

    return MSI_calls


def CheckForStutter(peaks):


    for i in range(0,len(peaks)-2):

        #if the
        if ((peaks[0][i+1]-peaks[0][i]) == 2) \
            and ((peaks[0][i+2]-peaks[0][i+1]) == 2):

            print("stutter detected")
            #thats possible stutter

            # now check if increasing or descreasing..
            #if ((peaks[0][i + 1] - peaks[0][i]) == 2) \
            #        and ((peaks[0][i + 2] - peaks[0][i + 1]) == 2):

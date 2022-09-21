import numpy as np
from scipy.stats import mode
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
import visuals

Liz500 = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500]

def getLadderPeaks(runName, trace_data_dictionary):

    #'DATA105' is the ladder channel
    ladder_trace = trace_data_dictionary['DATA105']

    highest_peaks_tup, smoothed_trace, threshold = findTop30PeaksLargestFirst(ladder_trace)

    #of the remaining peaks, keep the greatest.
    highest_tup=highest_peaks_tup[0]

    #now, keep the best 15 startibng from the right-most index.
    highest_peaks_tup.sort(key=lambda x: x[0], reverse=True)
    right_most = highest_peaks_tup[0:15]

    #put it together
    sixteen_peaks=[highest_tup] + right_most

    #want your ladder peaks leftmost on the left! not sorted by size
    sixteen_peaks.sort(key=lambda x: x[0])

    #bake LIZ500 into the peak tuple, for use later on.
    sixteen_peaks= [ (sixteen_peaks[i][0],sixteen_peaks[i][1],Liz500[i])
                       for i in range(0,16)]

    visuals.plotLadder(runName, threshold, smoothed_trace,
                       sixteen_peaks)

    return sixteen_peaks, threshold


def findTop30PeaksLargestFirst(ladder_trace):

    # very basic smoothing
    kernel_size = 20
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_trace = np.convolve(ladder_trace, kernel, mode='same')

    # calculate the threshold:
    ladder_variance = np.var(smoothed_trace)
    ladder_mode = mode(smoothed_trace)[0][0]
    ladder_sigma = np.sqrt(ladder_variance)
    threshold = ladder_mode + 0.25 * ladder_sigma

    # https://plotly.com/python/peak-finding/
    indices = find_peaks(smoothed_trace, height=threshold, distance=10)[0]
    peak_heights_tup = [(x, smoothed_trace[x]) for x in indices]
    # sort, greatest peak height first. Keep the best 30
    peak_heights_tup.sort(key=lambda x: x[1], reverse=True)
    highest_peaks_tup = peak_heights_tup[0:30]

    return highest_peaks_tup, smoothed_trace, threshold


def buildInterpolationdBasedOnLadder(run_folder, sixteen_peaks):

    #you are building a mapping from A -> B.
    #A = the peak positions in raw gel-travellijng space
    #B = the peak positions in bp length, as determined by the ladder.

    A = [x[0] for x in sixteen_peaks] #get the peak position, not intensity
    B = [x[2] for x in sixteen_peaks] #the LIZ 500 we baked in

    f = CubicSpline(A, B, bc_type='natural')

    #dont appy the spline where we dont have data!
    left_domain_limit = sixteen_peaks[0][0] - 100
    right_domain_limit = sixteen_peaks[15][0] + 100

    #test it looks good
    visuals.plotMapping(run_folder, A, B, left_domain_limit, right_domain_limit, f)

    return f, left_domain_limit, right_domain_limit

def remapATrace(tracedata_x_coords, tracedata_y_coords,
                fxn,  left_domain_limit, right_domain_limit):

    x_raw=[]
    y_raw=[]

    #get the subset of coords in the domain
    for i in range(0,len(tracedata_x_coords)):
        x = tracedata_x_coords[i]
        if  x >= left_domain_limit and x <= right_domain_limit:
            x_raw.append(tracedata_x_coords[i])
            y_raw.append(tracedata_y_coords[i])

    if fxn:
        x_new = fxn(x_raw)
    else:
        x_new= x_raw

    return x_raw, x_new, y_raw #y_raw and y_new are the same


def RemapLadder(run_folder, trace_data_dictionary, fxn,
                          left_domain_limit, right_domain_limit,
                          sixteen_peaks, threshold ):

    original_ladder=trace_data_dictionary['DATA105']
    num_data_points=len(original_ladder)
    old_x_coords=[x for x in range(0,num_data_points)]
    x_raw, x_new, y_new = remapATrace(old_x_coords,
                                      original_ladder,
                                      fxn, left_domain_limit, right_domain_limit)

    peak_xs=Liz500
    peak_ys= [peak[1] for peak in sixteen_peaks]

    visuals.plotRemappedTrace(run_folder,
                              x_new,
                              y_new , peak_xs, peak_ys,
                              threshold, "Ladder Trace", "Ladder Trace",
                              "Ladder Trace" , "Remapped")

    return x_new


def RemapDataTrace(run_folder, relevant_loci,
                          trace_data_dictionary, fxn,
                          left_domain_limit, right_domain_limit,
                          sixteen_peaks, channel_number ):

    channel_name= 'DATA' + str(channel_number)
    wavelength = trace_data_dictionary['DyeW' + str(channel_number)]
    channel_dye_name = trace_data_dictionary['DyeN' + str(channel_number)]
    raw_trace_data = (trace_data_dictionary[channel_name])


    highest_peaks_tup, smoothed_trace, threshold = findTop30PeaksLargestFirst(raw_trace_data)

    peak_xs = [peak[0] for peak in highest_peaks_tup]
    peak_ys = [peak[1] for peak in highest_peaks_tup]

    #plot raw trace, smoothed trace, and discovered peaks
    visuals.plotUnmappedTraceByColor(
        run_folder,
        raw_trace_data,smoothed_trace,highest_peaks_tup,
        wavelength, channel_dye_name, channel_number,"Raw")


    num_data_points=len(smoothed_trace)
    old_x_coords=[x for x in range(0,num_data_points)]
    trace_x_raw, trace_x_new, trace_y_new = remapATrace(old_x_coords, smoothed_trace,
                                      fxn, left_domain_limit, right_domain_limit)

    # if peaks are at the extreme end of the ladder, throw them out
    peak_x_raw, peak_x_new, peak_y_new = remapATrace(peak_xs,
                                      peak_ys,
                                      fxn,
                                      sixteen_peaks[1][0] + 10,
                                      sixteen_peaks[14][0] - 10)

    print("old peaks: " + str(peak_xs) )
    print("new peaks: " + str(peak_x_new) )
    plot_domain = [trace_x_new[0], trace_x_new[len(trace_x_new)-1]]

    visuals.plotRemappedTrace(run_folder, trace_x_new, trace_y_new,
                              peak_x_new, peak_y_new,
                              threshold,
                              wavelength, channel_dye_name, channel_number, "Remapped", plot_domain)

    # Plot a zoomed-in view for each loci of interest,
    # for this dye (usually one or two loci per dye).
    # Save off the MSI calls for each plot

    Peaks_inside_loci={}
    for loci_name in relevant_loci.keys():

        loci = relevant_loci[loci_name]
        dye_in_panel = loci["dye"]
        if dye_in_panel == channel_dye_name:
            print("Scanning " +  dye_in_panel + " trace")

        elif (dye_in_panel == "FAM" ) and \
                (channel_dye_name== "6-FAM"):
            #ok, good enough.
            print("Scanning FAM trace")

        else:
            continue

        bp_start = loci["length"][0] - 20
        bp_end = loci["length"][1] + 20
        plot_domain = [bp_start, bp_end ]

        visuals.plotRemappedTrace(run_folder, trace_x_new, trace_y_new,
                                  peak_x_new, peak_y_new,
                                  threshold,
                                  wavelength, channel_dye_name, channel_number,
                                  loci_name + "_Remapped",  plot_domain)

        peaks_in_loci_range=[]
        for i in range(0, len(peak_x_new)):
            x = peak_x_new[i]
            y = peak_y_new[i]
            if  ((x >=  loci["length"][0]) and (x <=  loci["length"][1])):
                peaks_in_loci_range.append([x,y])

        peaks_in_loci_range.sort(key=lambda x: x[0])
        Peaks_inside_loci[loci_name] = peaks_in_loci_range


    return Peaks_inside_loci, trace_x_new, trace_y_new

def PeaksToMsiCalls(peaks, trace_x_new, trace_y_new):

    if len(peaks) <= 1:
        consolidated_MSI=peaks
    else:
        #consoldate, if multiple peaks within one BP
        consolidated_MSI=[]
        i = 0
        while i < (len(peaks)-1):

            xi0=peaks[0][i]
            yi0=peaks[1][i]
            xi1=peaks[0][i+1]
            yi1=peaks[1][i+1]

            if (xi1 - xi0) < 1:
                consolidated_MSI.append( [ xi1 - xi0 / 2.0 , yi1 - yi0 / 2.0] )
                i = i +2
            else:
                consolidated_MSI.append(peaks[i])
                i = i +1
                if i == len(peaks)-1:
                    consolidated_MSI.append(peaks[i])

    #Check for stutter - do we have a sequence of peaks 2bp apart?
    #Is the height always increasing over the stutter?
    #if yes, skip the lower (left) peaks, and take the highest (right-most)



    MSI_calls = [round(peak) for peak in consolidated_MSI]
    MSI_calls = [ round(peak,1) for peak in consolidated_MSI]

    return MSI_calls
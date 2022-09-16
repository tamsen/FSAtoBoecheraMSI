import numpy as np
from scipy.stats import mode
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
import matplotlib
import matplotlib.pyplot as plt

import visuals

Liz500 = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500]

def getLadderPeaks(trace_data_dictionary):

    #'DATA105' is the ladder channel
    ladder_trace = trace_data_dictionary['DATA105']

    #very basic smoothing
    kernel_size = 20
    kernel = np.ones(kernel_size) / kernel_size
    smoothed_trace = np.convolve(ladder_trace , kernel, mode='same')

    #calculate the threshold:
    ladder_variance = np.var(smoothed_trace)
    ladder_mode = mode(smoothed_trace)[0][0]
    ladder_sigma = np.sqrt(ladder_variance)
    threshold = ladder_mode + 0.25*ladder_sigma

    # https://plotly.com/python/peak-finding/
    indices = find_peaks(smoothed_trace, height=threshold,distance=10)[0]
    peak_heights_tup= [ (x,smoothed_trace[x]) for x in indices ]

    #sort, greatest peak height first. Keep the best 30
    peak_heights_tup.sort(key=lambda x: x[1], reverse=True)
    highest_peaks_tup=peak_heights_tup[0:30]

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

    visuals.plotLadder(trace_data_dictionary, threshold, smoothed_trace,
                       sixteen_peaks)

    return sixteen_peaks

def buildInterpolationdBasedOnLadder(sixteen_peaks):

    #you are building a mapping from A -> B.
    #A = the peak positions in raw gel-travellijng space
    #B = the peak positions in bp length, as determined by the ladder.

    A = [x[0] for x in sixteen_peaks] #get the peak position, not intensity
    B = [x[2] for x in sixteen_peaks] #the LIZ 500 we baked in

    f = CubicSpline(A, B, bc_type='natural')

    print("A=" + str(A))
    print("B=" + str(B))

    fig, ax = plt.subplots(figsize=(10, 10))

    plt.plot(A, B, )
    # plt.plot(peak_indexes,peak_heights,"*")
    # plt.plot([threshold for x in new_x], "-", color="g")

    plt.title("MySpline")
    plt.xlabel("A = raw")
    plt.ylabel("B = remapped")

    # plt.show()
    # ax.legend(loc="upper right", title="Legend")
    plt.savefig("./tmp/mapping" + ".png")
    plt.close()


    return f

def remapATrace(tracedata_x_coords, fxn):

    new_x_coords = fxn(tracedata_x_coords)

    return new_x_coords


def testRemappingOnLadder(trace_data_dictionary, fxn, threshold):

    original_ladder=trace_data_dictionary['DATA105']
    num_data_points=len(original_ladder)
    old_x_coords=[x for x in range(0,num_data_points)]
    new_x_coords = remapATrace(old_x_coords, fxn)

    #cut off anything below zero
    x_to_plot=[]
    y_to_plot=[]

    for i in range(0,num_data_points):

        x=new_x_coords[i]
        y=original_ladder[i]
        if x >= 0 and x <= 100:
            x_to_plot.append(x)
            y_to_plot.append(y)


    #x datapoints start getting smaller again, so something seems wrong with the spline
    print("x=" + str(x_to_plot))
    print("y=" + str(y_to_plot))

    visuals.plotRemappedTrace(x_to_plot,
                              y_to_plot, threshold)

    return new_x_coords

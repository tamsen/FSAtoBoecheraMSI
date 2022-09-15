import numpy as np
from scipy.stats import mode
from scipy.signal import find_peaks
import visuals

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

    visuals.plotLadder(trace_data_dictionary, threshold, smoothed_trace,
                       sixteen_peaks)

    return sixteen_peaks
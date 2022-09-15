import matplotlib
import matplotlib.pyplot as plt

def plotLadder(trace_data_dictionary, threshold,
               smoothed_ladder_data, sixteen_peaks):

    peak_heights = [ x[1] for x in sixteen_peaks]
    peak_indexes = [x[0] for x in sixteen_peaks]


    fig, ax = plt.subplots(figsize=(10, 10))

    plt.plot(smoothed_ladder_data)
    plt.plot(peak_indexes,peak_heights,"*")
    plt.plot([threshold for x in smoothed_ladder_data], "-", color="g")

    first_peak=sixteen_peaks[0]
    ax.text(first_peak[0], first_peak[1], str(int(first_peak[0])))
    for peak in sixteen_peaks[1:16]:
        ax.text(peak[0], peak[1], str(int(peak[0])), rotation=45)


    plt.title("Ladder Trace")
    plt.xlabel("Distance Fragment Travelled")
    plt.ylabel("Intensity")

    sixteen_peaks.sort(key=lambda x: x[0])
    ladder_string=str([x[0] for x in sixteen_peaks])

    ax.text(20, 20000,
            "Ladder positions: \n\n" + ladder_string, style='italic',
            bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})

    #plt.show()
    #ax.legend(loc="upper right", title="Legend")
    plt.savefig("./tmp/ladder_plot" + ".png")
    plt.close()

def plotRawTraceByColor(trace_data_dictionary):


    wavelength = trace_data_dictionary['DyeW1']
    dyename = trace_data_dictionary['DyeN1']

    fig, ax = plt.subplots(figsize=(10, 10))

    text_overlay="wavelength: " + str(wavelength ) + "\n\ndye name: " \
                 + str(dyename)
    ax.text(6000, 20000, text_overlay, style='italic',
            bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})



    #'DATA1' is the first data channel
    plt.plot(trace_data_dictionary['DATA1'])

    plt.title("Dye1 Trace")
    plt.xlabel("Distance Travelled")
    plt.ylabel("Intensity")

    #plt.show()
    plt.savefig("./tmp/data1_plot" + ".png")
    plt.close()

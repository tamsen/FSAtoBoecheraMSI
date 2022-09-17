import numpy as np
import matplotlib.pyplot as plt

def plotLadder(run_folder, threshold,
               smoothed_ladder_data, sixteen_peaks):

    peak_heights = [ x[1] for x in sixteen_peaks]
    peak_indexes = [x[0] for x in sixteen_peaks]


    fig, ax = plt.subplots(figsize=(10, 10))

    plt.plot(smoothed_ladder_data)
    plt.plot(peak_indexes,peak_heights,"*")
    plt.plot([threshold for x in smoothed_ladder_data], "-", color="g")

    for peak in sixteen_peaks:

        ax.text(peak[0], -2000, "raw="+str(int(peak[0])), rotation=-45)
        ax.text(peak[0], peak[1]+500, "BP="+str(peak[2]), rotation=45)

    plt.title("Ladder Trace")
    plt.xlabel("Distance Fragment Travelled (raw)")
    plt.ylabel("Intensity")

    sixteen_peaks.sort(key=lambda x: x[0])
    ladder_string=str([x[0] for x in sixteen_peaks])
    bp_string=str([x[2] for x in sixteen_peaks])

    ax.text(20, 20000,
            "Ladder positions: \n\n" + ladder_string, style='italic',
            bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})

    ax.text(20, 15000,
            "BP lengths: \n\n" + bp_string, style='italic',
            bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})

    plt.savefig(run_folder  + "/Raw_Ladder_plot" + ".png")
    plt.close()


def plotRemappedTrace(run_folder, new_x,old_y, peak_xs,  peak_ys, threshold,
                      wavelength, dyename, channel_number, plot_prefix):

    fig, ax = plt.subplots(figsize=(10, 10))

    plt.plot(new_x,old_y)
    plt.plot(new_x,[threshold for x in new_x], "-", color="g")

    for i in range(0,len(peak_xs)):

         formatBPstring=  str('{0:3.1f}'.format(peak_xs[i]))
         ax.text(peak_xs[i], 6000, "BP="+str(formatBPstring), rotation=45)

    plt.title(plot_prefix + " " + dyename )
    plt.xlabel("Distance Fragment Travelled (BP)")
    plt.ylabel("Intensity")

    #ax.legend(loc="upper right", title="Legend")
    plt.savefig(run_folder + "/" + plot_prefix + "_" + dyename + "_plot" + ".png")

    plt.close()

def plotMapping(run_folder, A, B, left_domain_limit, right_domain_limit, f):

    num_data_points = (right_domain_limit - left_domain_limit) * .1 + 1.0
    x_original = list(np.linspace( left_domain_limit, right_domain_limit, int(num_data_points)))
    x_new = f(x_original)

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(A, B, "*", color="brown", label='Ladder Points')
    plt.plot(x_original, x_new, "-", color="blue", label='Cubic Spline Interpolation')
    plt.title("Mapping According To Ladder")
    plt.xlabel("A = raw distance travelled in gel")
    plt.ylabel("B = basepair equivalent")
    ax.legend(loc="upper left", title="Legend")
    plt.savefig(run_folder + "/Raw_to_BP_mapping" + ".png")
    plt.close()

def plotUnmappedTraceByColor(run_folder, trace_data,smoothed_trace,highest_peaks_tup,
                             wavelength, dyename, channel_number, plot_prefix):

    fig, ax = plt.subplots(figsize=(10, 10))

    text_overlay="wavelength: " + str(wavelength ) + "\n\ndye name: " \
                 + str(dyename)
    ax.text(6000, 20000, text_overlay, style='italic',
            bbox={'facecolor': 'gray', 'alpha': 0.5, 'pad': 10})

    peak_xs=[peak[0] for peak in highest_peaks_tup]
    peak_ys=[peak[1] for peak in highest_peaks_tup]

    #'DATA1' is the first data channel
    plt.plot(trace_data,  label=plot_prefix + ' Data')
    plt.plot(smoothed_trace, "-" , label='Smoothed Data')

    plt.plot( peak_xs, peak_ys, "*", label='Peaks')

    plt.title( plot_prefix + " " + dyename + " Trace")
    plt.xlabel("Distance Travelled")
    plt.ylabel("Intensity")

    ax.legend(loc="upper right", title="Legend")
    plt.savefig(run_folder +"/" + plot_prefix + "_" + dyename + "_plot" + ".png")
    plt.close()

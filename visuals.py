import numpy as np
import matplotlib.pyplot as plt

def plot_ladder(run_folder, threshold,
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


def plot_remapped_trace(run_folder, new_x, old_y, peak_xs, peak_ys, threshold,
                        wavelength, dyename, channel_number,
                        plot_prefix, plot_domain=False):

    #fig, ax = plt.subplots(figsize=(10, 10))
    fig, ax = plt.subplots()

    plt.plot(new_x, old_y)

    for i in range(0, len(peak_xs)):
        formatBPstring = str('{0:3.1f}'.format(peak_xs[i]))
        ax.text(peak_xs[i], peak_ys[i] + 100, "BP=" + str(formatBPstring), rotation=45)

    plt.plot(new_x, [threshold for x in new_x], "-", color="g")
    plt.title(plot_prefix + " " + dyename )
    plt.xlabel("Distance Fragment Travelled (BP)")
    plt.ylabel("Intensity")

    if (plot_domain):
        plt.xlim(plot_domain)

    if (len(peak_ys) > 0):
        y_max = max(peak_ys)
        plt.ylim([0,y_max+1000])

    #ax.legend(loc="upper right", title="Legend")
    plt.savefig(run_folder + "/" + plot_prefix + "_" + dyename + "_plot" + ".png")

    plt.close()


def plot_trace_and_special_points(fig, plot_index, new_x, new_y,
                                  peak_xs, peak_ys, loci, domain,dye_color, add_special_point_text):

    ax = fig.add_subplot(*plot_index)
    ax = plt.plot(new_x, new_y, c=dye_color)
    ax = plt.scatter(peak_xs, peak_ys, marker="*", c='orange')

    ax = plt.gca()
    ax = plt.text(0.05, 0.95, loci, horizontalalignment='left',
         verticalalignment='top', transform=ax.transAxes)

    if (add_special_point_text):
        ax = plt.gca()
        ax = plt.text(0.95, 0.95, "obs: " + str(peak_xs), horizontalalignment='right',
            verticalalignment='top', transform=ax.transAxes, fontsize=8)

    if (domain):
        plt.xlim(domain)

    if (len(peak_ys) > 0):
        y_max = max(peak_ys)
        plt.ylim([0,(y_max*1.10)])

    return ax

def write_per_sample_summary_plots(run_folder, by_sample_results, panel_info):
    
    # specific order to arrange the plots
    ordered_loci_list = ["dummy_index", "ICE3", "BF20", "A1",
                         "BF11", "ICE14", "C8",
                         "BF9", "BF18", "E9",
                         "BF3", "BF19", "B6",
                         "BF15", "Bdru266", "A3"]

    for sample_name, sample_result in by_sample_results.items():
      plot_traces_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)
      plot_ladders_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)
      plot_mappings_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)

def plot_ladders_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):
    
    fig = plt.figure(figsize=(10, 10))
    #ladder_plot_data = [runName, threshold, smoothed_trace, sixteen_peaks]

    for i in range(0,5):

        PS=i+1
        loci=ordered_loci_list[(i*3)+1]
        [runName, threshold, smoothed_trace, sixteen_peaks] = sample_result[loci].ladder_plotting_data
        domain=[500,sixteen_peaks[-1][0]+200]
        ax = plot_trace_and_special_points(fig, (5,1,PS),
                                           range(0, len(smoothed_trace)),
                                           smoothed_trace,
                                           [x[0] for x in sixteen_peaks],
                                           [x[1] for x in sixteen_peaks],
                                           "", domain, 'black', False)

        plt.ylabel("PS" + str(PS),  fontsize=16)

    fig.suptitle(sample_name + " all loci " )
    plt.savefig(run_folder + "/" + sample_name +"_ladder_for_loci"+ "_plot" + ".png")

    plt.close()

def plot_mappings_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):
    
    
    #mapping_plot_data = [run_folder, x_new, y_new, peak_xs, peak_ys, threshold, "Ladder Trace", "Ladder Trace",
    #                 "Ladder Trace", "Remapped"]
    
    fig = plt.figure(figsize=(10, 10))

    for i in range(1,16):

        loci=ordered_loci_list[i]
        [mapping_plot_data_spline, mapping_plot_data_linear] = sample_result[loci].mapping_plotting_data
        [x_original, x_new, A, B] = mapping_plot_data_spline

        ax = plot_trace_and_special_points(fig, (5,3,i), x_original, x_new, A, B,
                                           loci, False, "purple",  False)

        if i in [1,4,7,10,13,16]:
            plt.ylabel("PS" + str(int((i + 2.0 )/3.0)),  fontsize=16)

    fig.suptitle(sample_name + " all loci " )
    plt.savefig(run_folder + "/" + sample_name +"_mappings_for_loci"+ "_plot" + ".png")

    plt.close()

def plot_traces_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):

    dye_to_color={"FAM":"blue", "VIC": "green"}

    fig = plt.figure(figsize=(10, 10))

    for i in range(1,16):

        loci=ordered_loci_list[i]
        [new_x, new_y, peak_xs, peak_ys, plot_prefix, domain, dye] = sample_result[loci].plotting_data_evidence
        ax = plot_trace_and_special_points(fig, (5,3,i), new_x, new_y, peak_xs, peak_ys,
                                           loci, domain, dye_to_color[dye], True)

        if i in [1,4,7,10,13,16]:
            plt.ylabel("PS" + str(int((i + 2.0 )/3.0)),  fontsize=16)

    fig.suptitle(sample_name + " all loci " )
    plt.savefig(run_folder + "/" + sample_name +"_all_loci"+ "_plot" + ".png")

    plt.close()

def plotMapping(run_folder, A, B, left_domain_limit, right_domain_limit, f,
                interpolation_type):

    num_data_points = (right_domain_limit - left_domain_limit) * .1 + 1.0
    x_original = list(np.linspace( left_domain_limit, right_domain_limit, int(num_data_points)))
    x_new = f(x_original)

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(A, B, "*", color="brown", label='Ladder Points')
    plt.plot(x_original, x_new, "-", color="blue", label='Cubic Spline Interpolation')
    plt.title("Mapping According To Ladder (" + interpolation_type +  ")")
    plt.xlabel("A = raw distance travelled in gel")
    plt.ylabel("B = basepair equivalent")
    ax.legend(loc="upper left", title="Legend")
    plt.savefig(run_folder + "/Raw_to_BP_mapping_" + interpolation_type + ".png")
    plt.close()

    plotting_data=[x_original, x_new, A,B]
    #check mapping gives monotonically increasing results. Otherwise, this is *sus*
    monotonic_test = np.all((np.diff(x_new)) >=0)

    return monotonic_test,plotting_data

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


import matplotlib.pyplot as plt
import accuracy

def plot_trace_and_special_points(fig, plot_index, new_x, new_y,
                                  peak_xs, peak_ys, loci, domain, dye_color,
                                  msi_call_text, ladder_peak_txt):
    ax = fig.add_subplot(*plot_index)
    ax = plt.plot(new_x, new_y, c=dye_color)
    ax = plt.scatter(peak_xs, peak_ys, marker="*", c='orange')


    ax = plt.gca()
    ax = plt.text(0.05, 0.95, loci, horizontalalignment='left',
                  verticalalignment='top', transform=ax.transAxes)

    if (msi_call_text):
        ax = plt.gca()
        ax = plt.text(0.95, 0.95, "obs: " + str(peak_xs), horizontalalignment='right',
                      verticalalignment='top', transform=ax.transAxes, fontsize=8)

    #if (ladder_peak_txt != False):
    #    ax = plt.gca()
    #    for x in peak_xs:
    #        ax = plt.text(x 0.95, str(x), fontsize=8)

    if (domain):
        plt.xlim(domain)

    if (len(peak_ys) > 0):
        y_max = max(peak_ys)
        plt.ylim([0, (y_max * 1.10)])

    return ax


def write_per_sample_summary_plots(run_folder, by_sample_results, panel_info,
                                   truth_info):
    # specific order to arrange the plots
    ordered_loci_list = ["dummy_index", "ICE3", "BF20", "A1",
                         "BF11", "ICE14", "C8",
                         "BF9", "BF18", "E9",
                         "BF3", "BF19", "B6",
                         "BF15", "Bdru266", "A3"]

    for sample_name, sample_result in by_sample_results.items():
        truth_for_this_sample = accuracy.find_truth_for_this_sample(sample_name, truth_info)

        plot_traces_for_the_sample(run_folder, sample_name,
                                   sample_result, ordered_loci_list, truth_for_this_sample)

        plot_ladders_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)
        plot_mappings_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)


def plot_ladders_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):
    fig = plt.figure(figsize=(10, 10))

    for i in range(0, 5):
        PS = i + 1
        loci = ordered_loci_list[(i * 3) + 1]
        [runName, plot_name, threshold, smoothed_trace, sixteen_peaks] = sample_result[loci].ladder_plotting_data
        domain = [500, sixteen_peaks[-1][0] + 200]
        x_values = range(0, len(smoothed_trace))
        peak_xs=[x[0] for x in sixteen_peaks]
        ax = plot_trace_and_special_points(fig, (5, 1, PS),
                                           x_values,
                                           smoothed_trace,
                                           peak_xs,
                                           [x[1] for x in sixteen_peaks],
                                           "", domain, 'black', False, True)

        ax = plt.plot(x_values, [threshold for x in x_values], c="green")

        #ax = plt.gca()
        #for peak in peak_xs:
        #    ax.text(str(int(peak)), peak, -250, rotation=0)

        plt.ylabel("PS" + str(PS), fontsize=16)

    fig.suptitle(str(sample_name) + " all loci ")
    plt.savefig(str(run_folder) + "/" + str(sample_name) + "_ladder_for_loci" + "_plot" + ".png")

    plt.close()


def plot_mappings_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):

    fig = plt.figure(figsize=(10, 10))

    for i in range(1, 16):

        loci = ordered_loci_list[i]
        [mapping_plot_data_spline, mapping_plot_data_linear] = sample_result[loci].mapping_plotting_data
        [x_original, x_new, A, B] = mapping_plot_data_spline

        ax = plot_trace_and_special_points(fig, (5, 3, i), x_original, x_new, A, B,
                                           loci, False, "purple", False, False)

        if i in [1, 4, 7, 10, 13, 16]:
            plt.ylabel("PS" + str(int((i + 2.0) / 3.0)), fontsize=16)

    fig.suptitle(sample_name + " all loci ")
    plt.savefig(run_folder + "/" + sample_name + "_mappings_for_loci" + "_plot" + ".png")

    plt.close()


def plot_traces_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list,
                               truth_for_this_sample):
    dye_to_color = {"FAM": "blue", "VIC": "green"}

    fig = plt.figure(figsize=(10, 10))

    for i in range(1, 16):

        loci = ordered_loci_list[i]
        [new_x, new_y, peak_xs, peak_ys, plot_prefix, domain, dye] = sample_result[loci].plotting_data_evidence
        ax = plot_trace_and_special_points(fig, (5, 3, i), new_x, new_y, peak_xs, peak_ys,
                                           loci, domain, dye_to_color[dye], True, False)

        if truth_for_this_sample != False:
            truth = truth_for_this_sample[loci]
            ax = plt.gca()
            ax = plt.text(0.95, 0.85, "exp: " + str(truth), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=8)

        if i in [1, 4, 7, 10, 13, 16]:
            plt.ylabel("PS" + str(int((i + 2.0) / 3.0)), fontsize=16)

    fig.suptitle(sample_name + " all loci ")
    plt.savefig(run_folder + "/" + sample_name + "_all_loci" + "_plot" + ".png")

    plt.close()

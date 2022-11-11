import os
from signal_processing.ladder_analysis import GLOBAL_Liz500
import matplotlib.pyplot as plt



def plot_trace_and_special_points(fig, plot_index, xs, ys,
                                  peak_xs, peak_ys, loci, plot_domain, dye_color,
                                  msi_call_text, ladder_peak_txt):
    ax = fig.add_subplot(*plot_index)
    ax = plt.plot(xs, ys, c=dye_color)
    ax = plt.scatter(peak_xs, peak_ys, marker="*", c='orange')

    ax = plt.gca()
    ax = plt.text(0.05, 0.95, loci, horizontalalignment='left',
                  verticalalignment='top', transform=ax.transAxes)

    if (ladder_peak_txt):
        ax = plt.gca()
        for i in range(0, len(peak_xs)):
            ax = plt.text(peak_xs[i], peak_ys[i]*1.05, str(GLOBAL_Liz500[i]), rotation=0)

    ys_within_domain = []
    if (plot_domain):
        plt.xlim(plot_domain)
        ys_within_domain = [ys[i] for i in range(0, len(xs))
                            if plot_domain[0] < xs[i] < plot_domain[1]]
    if (len(peak_ys) > 0):
        y_max = max(peak_ys)
        plt.ylim([0, y_max * 1.1])

    if not ladder_peak_txt:
        if len(ys_within_domain) > 0:
            y_max = max(ys_within_domain)
            plt.ylim([0, y_max * 1.5])

    return ax


def write_per_sample_summary_plots(run_folder, by_sample_results):
    # specific order to arrange the plots
    ordered_loci_list = ["dummy_index", "ICE3", "BF20", "A1",
                         "BF11", "ICE14", "C8",
                         "BF9", "BF18", "E9",
                         "BF3", "BF19", "B6",
                         "BF15", "Bdru266", "A3"]

    for sample_name, sample_result in by_sample_results.items():
        # truth_for_this_sample = accuracy.find_truth_for_this_sample(sample_name, truth_info)

        plot_traces_for_the_sample(run_folder, sample_name,
                                   sample_result, ordered_loci_list)

        plot_ladders_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)
        plot_mappings_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list)


def plot_ladders_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):
    fig = plt.figure(figsize=(10, 10))

    for i in range(0, 5):
        PS = i + 1
        loci = ordered_loci_list[(i * 3) + 1]
        warning = ""
        if loci in sample_result.keys():
            [runName, plot_name, threshold, smoothed_trace, sixteen_peaks] = sample_result[loci].ladder_plotting_data
        else:
            # some dummy todata
            threshold = 0
            smoothed_trace = [1, 2, 3]
            sixteen_peaks = [[1, 2], [2, 3]]
            warning = "this loci had no ladder data"

        domain = [500, sixteen_peaks[-1][0] + 200]
        x_values = range(0, len(smoothed_trace))
        peak_xs = [x[0] for x in sixteen_peaks]
        ax = plot_trace_and_special_points(fig, (5, 1, PS),
                                           x_values,
                                           smoothed_trace,
                                           peak_xs,
                                           [x[1] for x in sixteen_peaks],
                                           warning, domain, 'black', False, True)

        ax = plt.plot(x_values, [threshold for x in x_values], c="green")

        plt.ylabel("PS" + str(PS), fontsize=16)

    fig.suptitle(str(sample_name) + " ladder (x=fragment travel distance, y=intensity) ")

    out_path = os.path.join(run_folder, "ladders")
    if not (os.path.exists(out_path)):
        os.makedirs(out_path)
    plt.savefig(out_path + "/" + str(sample_name) + "_ladder_for_loci" + "_plot" + ".png")

    plt.close()


def plot_mappings_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):
    fig = plt.figure(figsize=(10, 10))

    for i in range(1, 16):

        loci = ordered_loci_list[i]

        if loci in sample_result.keys():
            [mapping_plot_data_spline, mapping_plot_data_linear] = sample_result[loci].mapping_plotting_data

            if mapping_plot_data_spline:
                [x_original, x_new, A, B] = mapping_plot_data_spline

                ax = plot_trace_and_special_points(fig, (5, 3, i), x_original, x_new, A, B,
                                                   loci, False, "purple", False, False)

            elif mapping_plot_data_linear:
                [x_original, x_new, A, B] = mapping_plot_data_linear

                ax = plot_trace_and_special_points(fig, (5, 3, i), x_original, x_new, A, B,
                                               loci, False, "purple", False, False)
            else:
                ax = plot_trace_and_special_points(fig, (5, 3, i), [1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3],
                                                   loci + "Failed", False, "purple", False, False)

        else:
            ax = plot_trace_and_special_points(fig, (5, 3, i), [1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3],
                                               loci + "Failed", False, "purple", False, False)

        if i in [1, 4, 7, 10, 13, 16]:
            plt.ylabel("PS" + str(int((i + 2.0) / 3.0)), fontsize=16)

    fig.suptitle(sample_name + " mapping for all loci (x=base pair length estimate, y= fragment travel distance)")

    out_path = os.path.join(run_folder, "mappings")
    if not (os.path.exists(out_path)):
        os.makedirs(out_path)
    plt.savefig(out_path + "/" + sample_name + "_mappings_for_loci" + "_plot" + ".png")

    plt.close()


def plot_traces_for_the_sample(run_folder, sample_name, sample_result, ordered_loci_list):
    dye_to_color = {"FAM": "blue", "VIC": "green"}
    fig = plt.figure(figsize=(10, 10))
    species = False

    for i in range(1, 16):

        loci = ordered_loci_list[i]

        if loci not in sample_result.keys():
            dye = "FAM"
            ax = plot_trace_and_special_points(fig, (5, 3, i), [1], [1], [1], [1],
                                               loci, [0, 1], dye_to_color[dye], True, False)
        else:
            [new_x, new_y, raw_peak_xs, raw_peak_ys, filtered_peak_xs, filtered_peak_ys,
             threshold, plot_prefix, domain, dye] = \
                sample_result[loci].plotting_data_evidence

            ax = plot_trace_and_special_points(fig, (5, 3, i), new_x, new_y, raw_peak_xs, raw_peak_ys,
                                               loci, domain, dye_to_color[dye], True, False)

            ax = plt.scatter(filtered_peak_xs, filtered_peak_ys, s=80, facecolors='none', edgecolors='r')
            ax = plt.plot(new_x, [threshold for x in new_x], c="gray")

        if loci in sample_result:
            loci_result=sample_result[loci]
            truth = loci_result.truth_data
            species = loci_result.true_species
            final_accuracy = loci_result.final_accuracy
            raw_accuracy = loci_result.raw_accuracy
            raw_alleles_called  = loci_result.raw_alleles_called
            filtered_alleles_called  = loci_result.filtered_alleles_called
            final_alleles_called = loci_result.final_alleles_called
        else:
            truth = False


        ax = plt.gca()
        base_text_height=0.95
        font_size=6
        if len(raw_alleles_called) < 6:
            ax = plt.text(0.95, 0.95, "raw: " + str(raw_alleles_called), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=font_size)
        else:
            base_text_height = 0.85
            ax = plt.text(0.95, 0.95, "raw: " + str(raw_alleles_called[0:4]),
                                  horizontalalignment='right',
                                  verticalalignment='top', transform=ax.transAxes, fontsize=font_size)
        ax = plt.gca()
        ax = plt.text(0.95, 0.85, "+ " + str(raw_alleles_called[4:-1]),
                                  horizontalalignment='right',
                                  verticalalignment='top', transform=ax.transAxes, fontsize=font_size)

        ax = plt.gca()
        ax = plt.text(0.95,base_text_height-0.1, "filtered: " + str(filtered_alleles_called), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=font_size)

        ax = plt.gca()
        ax = plt.text(0.95, base_text_height-0.2, "final: " + str(final_alleles_called ), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=font_size)
        if truth:
            ax = plt.gca()
            ax = plt.text(0.95, base_text_height-0.3, "exp: " + str(truth), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=font_size)
            ax = plt.gca()
            ax = plt.text(0.95, base_text_height-0.4, "raw acc.: " + str(raw_accuracy), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=font_size)

            ax = plt.gca()
            ax = plt.text(0.95,base_text_height-0.5, "final acc.: " + str(final_accuracy), horizontalalignment='right',
                          verticalalignment='top', transform=ax.transAxes, fontsize=font_size)


        if i in [1, 4, 7, 10, 13, 16]:
            plt.ylabel("PS" + str(int((i + 2.0) / 3.0)), fontsize=16)


    plot_title = sample_name + " all loci"
    if species:
        plot_title = plot_title + "\n" + species
    fig.suptitle(plot_title)

    out_path= os.path.join(run_folder,"allele_calls")
    if not (os.path.exists(out_path)):
        os.makedirs(out_path)
    plt.savefig(out_path + "/" + sample_name + "_all_loci" + "_plot" + ".png")

    plt.close()

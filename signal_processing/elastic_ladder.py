import statistics

from signal_processing.ladder_analysis import GLOBAL_Liz500


def build_elastic_ladder_from_right(highest_peaks_with_index_from_right, confident_right_peak, num_peaks_needed):

    tolerances = get_tolerances()
    index_of_last_ladder_peak = confident_right_peak[2]
    ladder_peaks = [confident_right_peak]

    for i in range(1, num_peaks_needed+1):  # 0-14. skip the one at the start of the ladder.
        print("i:" + str(i))
        ith_peak_to_find = num_peaks_needed + 1 - i
        print("ith peak to find:" + str(ith_peak_to_find))
        print("ladder peaks so far:" + str(ladder_peaks))
        next_peak = get_peak_i(highest_peaks_with_index_from_right, ith_peak_to_find,
                               index_of_last_ladder_peak, tolerances)
        index_of_last_ladder_peak = next_peak[2]
        ladder_peaks = [next_peak] + ladder_peaks

    print("final ladder peaks:" + str(ladder_peaks))
    return ladder_peaks

def build_elastic_ladder(highest_peaks_tup):
    # get the right indexes
    highest_peaks_tup.sort(key=lambda x: x[0])
    highest_peaks_with_index = [[highest_peaks_tup[i][0], highest_peaks_tup[i][1], i] for i in
                                range(0, len(highest_peaks_tup))]
    tolerances = get_tolerances()
    # grab the right-most, big peak:
    # now, keep the best 16 starting from the right-most index.
    highest_peaks_with_index.sort(key=lambda x: x[0], reverse=True)
    sixteenth_peak = get_peak_500_(highest_peaks_with_index)
    sixteenth_peak_index = sixteenth_peak[2]

    highest_peaks_with_index.sort(key=lambda x: x[0])
    index_of_last_ladder_peak = sixteenth_peak_index
    ladder_peaks = [sixteenth_peak]
    num_peaks = 16
    # print("highest_peaks_tup:" + str(highest_peaks_tup))

    for i in range(1, num_peaks):  # 0-15
        print("i:" + str(i))
        ith_peak_to_find = num_peaks - i
        print("ith peak to find:" + str(ith_peak_to_find))
        print("ladder peaks so far:" + str(ladder_peaks))
        next_peak = get_peak_i(highest_peaks_with_index, ith_peak_to_find, index_of_last_ladder_peak, tolerances)
        index_of_last_ladder_peak = next_peak[2]
        ladder_peaks = [next_peak] + ladder_peaks

    print("final ladder peaks:" + str(ladder_peaks))
    return ladder_peaks


def get_tolerances():
    # labels=[50-35,75-50,100-75,139-100,150-139,160-150,200-160,250-200,300-250,340-300,350-340,400-350,
    #                450-400,490-450,500-490]

    # mean,range_min,range_max
    expected_x_diffs = [[145, 142.6, 147.4], [276.66, 272.3, 280.7], [265.5, 262.5, 268.5], [434.25, 426.8, 441.2],
                        [109.92, 108.8, 111.2], [111.33, 109.7, 113.3], [460.75, 453.7, 469.3], [563.92, 558.9, 572.1],
                        [629.92, 614.7, 642.3], [468.92, 461.8, 476.2], [122.67, 117.1, 127.9], [629.08, 611, 647],
                        [585.58, 566.8, 605.2], [488.92, 469.6, 510.4], [100.75, 96.2, 105.8]]

    # range_min,range_max
    expected_Y_diffs = [[-610.15, 1.25], [16.871, 155.339], [-10.475, 191.425], [-9.164, 253.3785],
                        [-9.55, 159.05], [-95.06, 94.17], [13.66, 220.54], [-30.3, 152.1],
                        [-177, 207], [-51.420, 150.42], [12.949, 104.641], [15.587, 473.783],
                        [0.74, 376.841], [16.87, 410.83], [-51.71, 49.21]]

    return [expected_x_diffs, expected_Y_diffs]


#def get_peak_500(highest_peaks_with_index):

def get_peak_500_(highest_peaks_with_index):
    # grab the right-most, big peak:
    # now, keep the best 20 starting from the right-most index.
    highest_peaks_with_index.sort(key=lambda x: x[0], reverse=True)
    right_peaks = highest_peaks_with_index[0:18]
    print("right peaks" + str(right_peaks))

    # We expect "500" to be between x=6400 and x=6800
    candidate_500s = [p for p in right_peaks if 6400 < p[0] < 6800]
    print("peaks in right spot" + str(candidate_500s))
    heights = [p[1] for p in candidate_500s]
    mean_height = statistics.mean(heights)
    tall_ones_only = [p for p in candidate_500s if p[1] >= mean_height]
    tall_ones_only.sort(key=lambda x: x[0])
    print("candidate_500s" + str(tall_ones_only))
    choice = tall_ones_only[-1]
    print("candidate_500 choice" + str(choice))

    return choice

def get_peak_i(peaks_from_left_to_right, nth_ladder_peak_needed,
               index_of_last_ladder_peak_in_trace_space, tolerances):
    last_peak_to_the_right_x = peaks_from_left_to_right[index_of_last_ladder_peak_in_trace_space][0]
    last_peak_to_the_right_y = peaks_from_left_to_right[index_of_last_ladder_peak_in_trace_space][1]
    print("last_peak_to_the_right_x:" + str(last_peak_to_the_right_x))
    print("last_peak_to_the_right_y:" + str(last_peak_to_the_right_y))

    expected_x = last_peak_to_the_right_x - tolerances[0][nth_ladder_peak_needed - 1][0]
    tolerances_range_x_for_i = tolerances[0][nth_ladder_peak_needed - 1][1:3]
    tolerances_range_y_for_i = tolerances[1][nth_ladder_peak_needed - 1][0:2]

    print("tolerances_range_x_for_i:" + str(tolerances_range_x_for_i))
    print("tolerances_range_y_for_i:" + str(tolerances_range_y_for_i))

    #i_start = max([index_of_last_ladder_peak_in_trace_space - 6, 0])
    #i_end = max(index_of_last_ladder_peak_in_trace_space +1, 7)
    i_start = max([index_of_last_ladder_peak_in_trace_space + 1, 0])
    i_end = min(index_of_last_ladder_peak_in_trace_space + 2, len(peaks_from_left_to_right)-1)
    print("candidate_peak_indexes:" + str([i_start, i_end]))
    print("Boo!! two!!")
    candidate_next_peaks = peaks_from_left_to_right[i_start:i_end]
    print("candidate_next_peaks:" + str(candidate_next_peaks))

    expected_x_range = [last_peak_to_the_right_x - tolerances_range_x_for_i[1],
                        last_peak_to_the_right_x - tolerances_range_x_for_i[0]]
    expected_y_range = [last_peak_to_the_right_y - tolerances_range_y_for_i[1],
                        last_peak_to_the_right_y - tolerances_range_y_for_i[0]]

    print("expected_x_range:" + str(expected_x_range))
    print("expected_y_range:" + str(expected_y_range))

    # how  many are within parameters?
    peaks_in_the_right_place = [p for p in candidate_next_peaks if
                                ((expected_x_range[0] < p[0] < expected_x_range[1]) and
                                 (expected_y_range[0] < p[1] < expected_y_range[1]))]

    print("peaks in right place:" + str(peaks_in_the_right_place))

    # if we found nothing, loosen up the height requirement
    if len(peaks_in_the_right_place) < 1:
        peaks_in_the_right_place = [p for p in candidate_next_peaks if
                                    (expected_x_range[0] < p[0] < expected_x_range[1])]

    # if we found nothing, throw the kitchen sink at it
    if len(peaks_in_the_right_place) < 1:
        peaks_in_the_right_place = candidate_next_peaks

    # All else fails, fake it till you make it?
    if len(peaks_in_the_right_place) == 0:
        return False
        # return[expected_x,-1]

    # no choice, thats it
    if len(peaks_in_the_right_place) == 1:
        return peaks_in_the_right_place[0]

    else:  # its good to have choices
        return get_peak_closest_to_x(expected_x, peaks_in_the_right_place)

    # All else fails, fake it till you make it?
    return False


def get_peak_closest_to_x(expected_x, peaks_in_the_right_place):
    # get the peak in the most likely spot, regardless of height
    print("expected_x:" + str(expected_x))
    observed_x = [peaks_in_the_right_place[i][0]
                  for i in range(0, len(peaks_in_the_right_place))]
    diffs_between_obs_and_exp_pos = [abs(observed_x[i] - expected_x)
                                     for i in range(0, len(peaks_in_the_right_place))]
    closest_it_comes = min(diffs_between_obs_and_exp_pos)
    for i in range(0, len(peaks_in_the_right_place)):
        d = abs(peaks_in_the_right_place[i][0] - expected_x)
        if d == closest_it_comes:
            return peaks_in_the_right_place[i]

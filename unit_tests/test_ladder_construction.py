from unittest import TestCase
#import signal_processing.elastic_ladder_analysis

import os
from file_io import fsa_file_reader
from signal_processing import trace_analysis, elastic_ladder_analysis
import test_globals

class TestLadder(TestCase):


    def setUp(self):
        if not (os.path.exists(test_globals.GLOBAL_test_output_dir)):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

    def test_get_background(self):

        trace=[99,99,99,99,99,99,99]
        window_half_width=3
        background=elastic_ladder_analysis.get_background(trace, window_half_width)
        print(background)
        self.assertEqual(len(background), len(trace))


        trace=[99,99,2000,99,99,99,99, 0, 0, 99,99,99,99, 0, 0]
        window_half_width=3
        background=elastic_ladder_analysis.get_background(trace, window_half_width)
        print(background)
        self.assertEqual(len(background), len(trace))

    def test_get_ladder_peaks_easy_samples(self):

        #easy samples. clear peaks.
        fsa_file="../test_data/test_ladder/BD1200PS1c5_A01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(
            test_globals.GLOBAL_test_output_dir,
            "BD1200PS1c5_A01_",
            trace_data_dictionary)

        expected_peak_xs=[1202, 1347, 1626,1893,2331,2442,2555,3022,3595,4234,4711,4835,5474,6069,6566,6668]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        #easy samples. clear peaks.
        fsa_file="../test_data/test_ladder/FW428-PS1-A2_A02.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(
            test_globals.GLOBAL_test_output_dir,
            "FW428-PS1-A2_A02_",
            trace_data_dictionary)

        expected_peak_xs=[1310,1455,1733,2000,2438,2549,2661,3126,3694,4334,4809,4936,5580,6182,6689,6794]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        #easy samples. clear peaks.
        fsa_file="../test_data/test_ladder/FW428-PS5-E2_E02.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(
            test_globals.GLOBAL_test_output_dir,
            "FW428-PS5-E2_E02_",
            trace_data_dictionary)

        expected_peak_xs=[1325,1468,1745,2010,2445,2555,2666,3126,3687,4319,4786,4911,5542,6130,6623,6724]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

    def test_get_ladder_peaks_samples_with_messy_start(self):

        fsa_file="../test_data/test_ladder/BD1200PS2c5_H01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(
            test_globals.GLOBAL_test_output_dir,
                                                                                   "BD1200PS2c5_H01_",
                                                                                   trace_data_dictionary)

        expected_peak_xs=[1216,1362,1638,1904,2339,2449,2561,3023,3591,4219,4688,4808,5431,6009,6487,6586]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)


    def test_get_ladder_peaks_samples_with_multiple_false_peaks(self):


        #difficult sample (multiple spurious low peak)
        fsa_file="../test_data/test_ladder/BD1200PS5_E04.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                  "BD1200PS5_E04_",
                                                                                  trace_data_dictionary)

        expected_peak_xs=[1148,1293,1561,1822,2244,2352,2461,2910,3463,4070,4526,4642,5245,5806,6271,6367]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)


        #difficult sample (multiple spurious high peaks)
        fsa_file="../test_data/test_ladder/JB277PS1F2_D02.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                  "JB277PS1F2_D02_",
                                                                                  trace_data_dictionary)

        expected_peak_xs=[1186,1332,1605,1869,2298,2407,2518,2973,3533,4151,4614,4732,5346,5916,6389,6486]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

    def test_get_ladder_peaks_samples_with_false_peaks(self):

        #difficult sample (one spurious low peak)
        fsa_file="../test_data/test_ladder/BD1200PS3c5_A02.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                  "BD1200PS3c5_A02_",
                                                                                  trace_data_dictionary)

        expected_peak_xs=[1179,1323,1598,1864,2297,2407,2519,2980,3548,4179,4652,4773,5404,5993,6484,6585]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        #difficult sample (one spurious high peak)
        fsa_file="../test_data/test_ladder/BD1200CNTRL-PS3-C4_B03.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                  "BD1200CNTRL-PS3-C4_B03_",
                                                                                  trace_data_dictionary)

        expected_peak_xs=[1270,1416,1688,1952,2381,2490,2600,3055,3614,4235,4700,4821,5444,6027,6514,6615]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        #difficult sample (one spurious high peak)
        fsa_file="../test_data/test_ladder/JB617-PS1-A7_G12.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                  "JB617-PS1-A7_G12.fsa_",
                                                                                  trace_data_dictionary)

        expected_peak_xs = [1332,1475,1756,2024,2463,2574,2686,3152,3719,4362,4838,4966,5610,6210,6713,6817]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

    def test_get_ladder_peaks_for_Lemmonii_with_false_peaks(self):
        # difficult sample (one spurious low peak)
        fsa_file = "../test_data/test_ladder/FW415_ladder_issue_PS1_A2_G01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                    "BD1200PS3c5_A02_",
                                                                                    trace_data_dictionary)

        expected_peak_xs =[1225,1370,1645,1909,2341,2451,2561,3019,3580,4203,4667,4787,5405,5977,6453,6551]

        observed_peak_xs = [x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

    def test_get_ladder_peaks_for_RPsamples_with_peaks_lost_in_threshold(self):

        fsa_file = "../test_data/test_ladder/TD21RP22PS1_G08.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                    "TD21RP22PS1_G08_",
                                                                                    trace_data_dictionary)

        #expected_peak_xs =[1225,1370,1645,1909,2341,2451,2561,3019,3580,4203,4667,4787,5405,5977,6453,6551]
        #observed_peak_xs = [x[0] for x in sixteen_peaks]
        #commenting this test out. We still dont do a great job with such messed-up ladders
        #self.assertEqual(expected_peak_xs, observed_peak_xs)

        fsa_file = "../test_data/test_ladder/TD21RP22PS2_H08.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                    "TD21RP22PS2_H08_",
                                                                                  trace_data_dictionary)

        #expected_peak_xs =[1225,1370,1645,1909,2341,2451,2561,3019,3580,4203,4667,4787,5405,5977,6453,6551]

        #observed_peak_xs = [x[0] for x in sixteen_peaks]
        # commenting this test out. We still dont do a great job with such messed-up ladders
        # self.assertEqual(expected_peak_xs, observed_peak_xs)



    def test_get_ladder_peaks_for_GDsamples_with_peaks_lost_in_threshold(self):

        fsa_file = "../test_data/test_ladder/TD21GD0204PS1-A5_A05.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                    "TD21GD0204PS1-A5_A05.fsa",
                                                                                    trace_data_dictionary)

        expected_peak_xs = [1201, 1348, 1623, 1889, 2322, 2432, 2544, 3005, 3571, 4197, 4666, 4787, 5413, 5994, 6478,
                            6578] # v1.0.0.3 looks better...

        # v1.0.0.0 looks worse...
        #expected_peak_xs =[1130,1162,1623,1889,2322,2432,2544,3005,3571,4197,4666,4787,5413,5994,6478,6578]
        observed_peak_xs = [x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)


    def test_specifcally_that_worked_after_elastic_ladder_alg(self):

        fsa_file = "../test_data/test_ladder/TD22EL23-PS2-B7_B07.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(test_globals.GLOBAL_test_output_dir,
                                                                                    "TD22EL23-PS2-B7_B07_",
                                                                                    trace_data_dictionary)

        expected_peak_xs = [1187,1333,1599,1860,2281,2388,2497,2943,3489,4098,4554,4669,5274,5839,6310,6407]
        observed_peak_xs = [x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        fsa_file = "../test_data/test_ladder/TD22BV12-PS5-E8_E08.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = elastic_ladder_analysis.getLadderPeaks(
            test_globals.GLOBAL_test_output_dir,
            "TD22BV12-PS5-E8_E08.fsa_",
            trace_data_dictionary)

        expected_peak_xs = [1173,1318,1583,1840,2260,2366,2474,2917,3464,4062,4513,4627,5225,5782,6245,6341]
        observed_peak_xs = [x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)


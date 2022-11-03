from unittest import TestCase
import signal_processing.ladder_analysis
import file_io.fsa_file_reader
import os

from file_io import fsa_file_reader
from signal_processing import trace_analysis, ladder_analysis

global test_output_dir


class Test(TestCase):


    def setUp(self):
        global test_output_dir
        test_output_dir="../test_case_output/"
        if not (os.path.exists(test_output_dir)):
            os.makedirs(test_output_dir)

    def test_get_ladder_peaks_easy_samples(self):

        #easy samples. clear peaks.
        fsa_file="../test_data/BD1200PS1c5_A01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = ladder_analysis.getLadderPeaks(test_output_dir,
                                                                                   "BD1200PS1c5_A01_",
                                                                                   trace_data_dictionary)

        expected_peak_xs=[1202,1347,1626,1893,2331,2442,2555,3022,3595,4234,4711,4835,5474,6069,6566,6668]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

    def test_get_ladder_peaks_samples_with_messy_start(self):

        fsa_file="../test_data/BD1200PS2c5_H01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = ladder_analysis.getLadderPeaks(test_output_dir,
                                                                                   "BD1200PS2c5_H01_",
                                                                                   trace_data_dictionary)

        expected_peak_xs=[1216,1362,1638,1904,2339, 2449,2561,3023,3591,4219,4688,4808,5431,6009,6487,6586]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)


    def test_get_ladder_peaks_samples_with_false_peaks(self):


        #difficult sample (one spurious low peak)
        fsa_file="../test_data/BD1200PS3c5_A02.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data =ladder_analysis.getLadderPeaks(test_output_dir,
                                                                                  "BD1200PS3c5_A02_",
                                                                                  trace_data_dictionary)
        expected_peak_xs=[ 1138,1324, 1599, 1864, 2299, 2408, 2521,
                           2982, 3550, 4178, 4653, 4775, 5405, 5993, 6484, 6586]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        #self.assertEqual(expected_peak_xs, observed_peak_xs)

        #difficult sample (one spurious high peak)
        fsa_file="../test_data/BD1200CNTRL-PS3-C4_B03.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data =ladder_analysis.getLadderPeaks(test_output_dir,
                                                                                  "BD1200CNTRL-PS3-C4_B03_",
                                                                                  trace_data_dictionary)

        expected_peak_xs=[1270,1416,1688,1952,2381,2490,2600,3055,3614,4235,4700,4821,5444,6027,6514,6615]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

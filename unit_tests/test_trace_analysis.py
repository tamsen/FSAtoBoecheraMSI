from unittest import TestCase
import signal_processing.trace_analysis
import file_io.fsa_file_reader
import os

from file_io import fsa_file_reader
from signal_processing import trace_analysis

global test_output_dir


class Test(TestCase):


    def setUp(self):
        global test_output_dir
        test_output_dir="../test_case_output/"
        if not (os.path.exists(test_output_dir)):
            os.makedirs(test_output_dir)

    def test_get_ladder_peaks(self):

        #easy samples. clear peaks.
        fsa_file="../test_data/BD1200PS1c5_A01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = trace_analysis.getLadderPeaks(test_output_dir,
                                                                                   "BD1200PS1c5_A01_",
                                                                                   trace_data_dictionary)
        expected_peak_xs=[1160, 1347, 1626, 1893, 2332, 2442, 2556, 3024, 3598, 4236, 4713, 4837, 5476, 6071, 6567, 6669]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        fsa_file="../test_data/BD1200PS2c5_H01.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data = trace_analysis.getLadderPeaks(test_output_dir,
                                                                                   "BD1200PS2c5_H01_",
                                                                                   trace_data_dictionary)
        expected_peak_xs=[1174, 1363,1639,1905, 2340,2450,2562, 3026,3594,4222,4691,4811,5434,6011,6490,6588]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)

        #difficult sample (one vey low peak)
        fsa_file="../test_data/BD1200PS3c5_A02.fsa"
        dye_to_channel_mapping, trace_data_dictionary = fsa_file_reader.readFSAFile(fsa_file)
        sixteen_peaks, threshold, ladder_plot_data =trace_analysis.getLadderPeaks(test_output_dir,
                                                                                  "BD1200PS3c5_A02_",
                                                                                  trace_data_dictionary)
        expected_peak_xs=[ 1324, 1599, 1864, 2065, 2299, 2408, 2521, 2982, 3550, 4178, 4653, 4775, 5405, 5993, 6484, 6586]
        observed_peak_xs=[x[0] for x in sixteen_peaks]
        self.assertEqual(expected_peak_xs, observed_peak_xs)



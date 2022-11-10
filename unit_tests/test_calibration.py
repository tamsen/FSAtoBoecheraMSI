import unittest
import os

import fsa_directory_processor
import test_globals
from file_io import xml_file_readers


class TestCalibration(unittest.TestCase):

    def setUp(self):
        if not (os.path.exists(test_globals.GLOBAL_test_output_dir)):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

    def test_LxPxR_calibration(self):

        path ="/home/tamsen/Data/Eton/LxPxR_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "calibrationTest",
                                                        "LxPxR_FSA_to_microsat_script_results")
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['BD1200'] >= 99.72, True)
        self.assertEqual(avg_sample_accuracy['BD1200CNTRL'] >= 99.00, True)
        self.assertEqual(avg_sample_accuracy['BP28'] >= 98.32, True) #brought down by BF9, PS3, calling extra peak.
        self.assertEqual(avg_sample_accuracy['LA846'] >= 99.44, True)


    def test_paupercula_calibration(self):

        path ="/home/tamsen/Data/Eton/Paupercula_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "calibrationTest",
                                                        "Paupercula_FSA_to_microsat_script_results")
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['FW346'] >= 88.98, True)
        self.assertEqual(avg_sample_accuracy['FW437'] >= 89.21, True)

    def test_lemmonii_calibration(self):

        path ="/home/tamsen/Data/Eton/Lemmonii_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "calibrationTest",
                                                        "Lemmonii_FSA_to_microsat_script_results")
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['FW428'] >= 82.33, True)
        self.assertEqual(avg_sample_accuracy['FW415'] >= 74.92, True)

    def test_retrofracta_calibration(self):

        path ="/home/tamsen/Data/Eton/Retrofracta_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "calibrationTest",
                                                        "Retrofracta_FSA_to_microsat_script_results")
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['JB276'] >= 75.05, True)
        self.assertEqual(avg_sample_accuracy['JB277'] >= 82.19, True)

if __name__ == '__main__':
    unittest.main()

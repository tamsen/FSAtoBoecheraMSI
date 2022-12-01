import unittest
import os
from datetime import datetime

import fsa_directory_processor
import test_globals
import version
from file_io import xml_file_readers


class TestCalibration(unittest.TestCase):


    def setUp(self):
        if not (os.path.exists(test_globals.GLOBAL_test_output_dir)):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

    def get_time_string(self):
        now = datetime.now()
        day = now.strftime("%d_%m_%Y")
        time = now.strftime("%H_%M_%S")
        time_stamp_string = "_".join([day, time])
        return time_stamp_string

    def test_LxPxR_calibration(self):
        path = "/home/tamsen/Data/Eton/LxPxR_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "CalibrationTest",
                                                        "LxPxR_FSA_to_microsat_script_results_" +
                                                        self.get_time_string())
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            version.version_info(),
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        #self.assertEqual(avg_sample_accuracy['BD1200'] >= 99.90, True)
        #self.assertEqual(avg_sample_accuracy['BD1200CNTRL'] >= 99.75, True)
        #self.assertEqual(avg_sample_accuracy['BP28'] >= 99.60, True) #fixed when we increase sutter-distance from 3.5 -> 4
        #self.assertEqual(avg_sample_accuracy['LA846'] >= 99.62, True)

        # Went down a bit when I noticed the BP28 truth was wrong.
        # I suspect there are errors in the original BP28 calls, which is why the BP28 score is low
        self.assertEqual(avg_sample_accuracy['BD1200'] >= 99.90, True)
        self.assertEqual(avg_sample_accuracy['BD1200CNTRL'] >= 99.90, True)
        self.assertEqual(avg_sample_accuracy['BP28'] >= 97.14, True) #fixed when we increase sutter-distance from 3.5 -> 4
        self.assertEqual(avg_sample_accuracy['LA846'] >= 99.75, True)

        # Also check the species designation are right
        BD1200_species_determination=by_sample_results['BD1200']['A1'].BMW_determination
        BD1200CNTRL_species_determination=by_sample_results['BD1200CNTRL']['A1'].BMW_determination
        BP28_species_determination=by_sample_results['BP28']['A1'].BMW_determination
        LA846_species_determination=by_sample_results['LA846']['A1'].BMW_determination

        self.assertEqual("lemmonii x paupercula x retrofracta" in BD1200_species_determination, True)
        self.assertEqual("lemmonii x paupercula x retrofracta" in BD1200CNTRL_species_determination, True)
        self.assertEqual("lemmonii x paupercula x retrofracta" in BP28_species_determination, True) #(not yet passing. TODO, fx this)
        self.assertEqual("lemmonii x paupercula x retrofracta" in LA846_species_determination, True)

    def test_paupercula_calibration(self):
        path = "/home/tamsen/Data/Eton/Paupercula_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "CalibrationTest",
                                                        "Paupercula_FSA_to_microsat_script_results_" +
                                                        self.get_time_string())
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            version.version_info(),
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['FW346'] >= 83.76, True)
        self.assertEqual(avg_sample_accuracy['FW437'] >= 94.49, True)
        #could be brough tup by fixing extra peak in BF3, but it breaks lxpxr accuracy
        #self.assertEqual(avg_sample_accuracy['FW437'] >= 91.84, True)

        FW346_species_determination=by_sample_results['FW346']['A1'].BMW_determination.split(" ")[0]
        FW437_species_determination=by_sample_results['FW437']['A1'].BMW_determination.split(" ")[0]

        self.assertEqual("paupercula" == FW346_species_determination, True)
        self.assertEqual("paupercula" == FW437_species_determination, True)

    def test_lemmonii_calibration(self):
        path = "/home/tamsen/Data/Eton/Lemmonii_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "CalibrationTest",
                                                        "Lemmonii_FSA_to_microsat_script_results_" +
                                                        self.get_time_string())
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        #by_sample_results[sample_name][loci] = msi_results_for_loci
        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            version.version_info(), all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['FW428'] >= 89.37, True)
        self.assertEqual(avg_sample_accuracy['FW415'] >= 85.01, True)

        FW428_species_determination=by_sample_results['FW428']['A1'].BMW_determination.split(" ")[0]
        FW415_species_determination=by_sample_results['FW415']['A1'].BMW_determination.split(" ")[0]

        self.assertEqual("lemmonii" == FW428_species_determination, True)
        self.assertEqual("lemmonii*" == FW415_species_determination, True)

    def test_retrofracta_calibration(self):
        path = "/home/tamsen/Data/Eton/Retrofracta_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "CalibrationTest",
                                                        "Retrofracta_FSA_to_microsat_script_results_" +
                                                        self.get_time_string())
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            version.version_info(),
            all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['JB276'] >= 77.14, True)
        self.assertEqual(avg_sample_accuracy['JB277'] >= 82.41, True)

        JB276_species_determination=by_sample_results['JB276']['A1'].BMW_determination.split(" ")[0]
        JB277_species_determination=by_sample_results['JB277']['A1'].BMW_determination.split(" ")[0]

        self.assertEqual("retrofracta" in JB276_species_determination, True)
        self.assertEqual("retrofracta" in JB277_species_determination, True)

    def test_miscellaneous_calibration(self):
        path = "/home/tamsen/Data/Eton/Miscellaneous_Calibration"
        truth_file = os.path.join("../data/TruthData.xml")
        panel_file = os.path.join("../data/Panel.xml")
        truth_info = xml_file_readers.read_truth_data(truth_file)
        panel_info = xml_file_readers.readPanelXml(panel_file)

        output_folder_inside_data_folder = os.path.join(test_globals.GLOBAL_test_output_dir,
                                                        "CalibrationTest",
                                                        "Miscellaneous_FSA_to_microsat_script_results_" +
                                                        self.get_time_string())
        results_specific_to_this_subfolder = {}
        all_results_by_file = {}

        by_sample_results, avg_loci_accuracy, avg_sample_accuracy = fsa_directory_processor.process_directory(
            version.version_info(), all_results_by_file, output_folder_inside_data_folder,
            panel_info, path, results_specific_to_this_subfolder, truth_info)

        self.assertEqual(avg_sample_accuracy['JB1488'] >= 78, True)  # f x f
        self.assertEqual(avg_sample_accuracy['JB617'] >= 70, True)  # l x n x r
        self.assertEqual(avg_sample_accuracy['FW1379'] >= 81, True)  # l x r
        self.assertEqual(avg_sample_accuracy['BP27'] >= 82, True)  # l x r


if __name__ == '__main__':
    unittest.main()

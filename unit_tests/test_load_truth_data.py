import unittest
import os
from file_io import xml_file_readers


class test_load_truth_data(unittest.TestCase):

    def test_load_calibration_truth(self):

        input_file_path= os.path.join("../data/TruthData.xml")
        truth_data_by_sample_name = xml_file_readers.read_truth_data(input_file_path)

        expected_num_calibration_samples = 13
        expected_num_loci_for_bd1200 = 15
        expected_species_name ='lxpxr'
        bd1200data=truth_data_by_sample_name["BD1200"]

        self.assertEqual(expected_num_calibration_samples, len(truth_data_by_sample_name.keys()))
        self.assertEqual(expected_num_loci_for_bd1200, len(bd1200data.sample_data_by_loci))
        self.assertEqual(expected_species_name, bd1200data.species_name)

if __name__ == '__main__':
    unittest.main()

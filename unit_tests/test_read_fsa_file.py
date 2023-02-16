import unittest
import os
from file_io import fsa_file_reader

class Test_Read_FSA_File(unittest.TestCase):

    def test_czech(self):

        fsa_file = os.path.join("../test_data","test_reading_czech_fsa_files", "B07736_PS3.fsa")
        dye_to_channel_mapping, all_collected_data = fsa_file_reader.readFSAFile(fsa_file)

        self.assertEqual(dye_to_channel_mapping['Dye1'], 1)
        self.assertEqual(all_collected_data['DATA1'][0], -1)
        self.assertEqual(all_collected_data['DATA4'][0], 11)


if __name__ == '__main__':
    unittest.main()

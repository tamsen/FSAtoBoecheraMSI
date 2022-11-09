import unittest
import os
import test_globals

class TestCalibration(unittest.TestCase):

    def setUp(self):
        if not (os.path.exists(test_globals.GLOBAL_test_output_dir)):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

    def test_calibration(self):
        #home / tamsen / Data / Eton / LxPxR_Calibration

        if os.path.exists(test_globals.GLOBAL_test_output_dir):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

        self.assertEqual(True, True,)


if __name__ == '__main__':
    unittest.main()

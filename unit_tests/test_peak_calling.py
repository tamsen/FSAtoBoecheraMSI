import unittest

from signal_processing import peak_analysis


class TestPeakFiltering(unittest.TestCase):

    def test_iteratively_clean_short_peaks(self):


        peaks = [[89.9, 823.1000000000001], [92.10233688946954, 855.5999999999999],\
               [93.9, 13956.1], [97.86468582458026, 1325.6], [99.90770391497513, 12466.0],\
               [101.83534326055506, 11292.900000000001], [103.01833113835109, 2365.2]]

        step_width_left = 10
        step_proportion_left = 0.4
        step_width_right = 3
        step_proportion_right = 0.4

        peaks = peak_analysis.iteratively_clean_short_peaks(peaks, step_width_left, step_proportion_left,
                                                            step_width_right, step_proportion_right)
        self.assertEqual(peaks[0][0], 93.9)


        peaks= [[102.6, 578.7], [121.4, 12071.3]]
        peaks = peak_analysis.iteratively_clean_short_peaks(peaks, 100, .1, 100, .1)
        self.assertEqual(peaks[0][0], 121.4)

        peaks= [[102.6, 578.7], [119.4, 1328.3], [121.4, 12071.3]]
        peaks = peak_analysis.iteratively_clean_short_peaks(peaks, 100, .1, 100, .1)
        self.assertEqual(len(peaks), 3)

        peaks= [[102.6, 578.7], [119.4, 1328.3], [121.4, 12071.3]]
        peaks = peak_analysis.iteratively_clean_short_peaks(peaks, 100, .2, 100, .2)
        self.assertEqual(len(peaks), 1)

    def test_brutally_clean_short_peaks(self):

        peaks = [[89.9, 1128.7], [95.6, 1314.4], [106.0, 1840.5], [108.0, 8654.2],
                 [109.8, 24494.600000000002], [110.9, 11339.9]]
        peaks = peak_analysis.brutally_clean_short_peaks(peaks, .2)
        print(str(peaks))
        self.assertEqual(len(peaks), 3)

if __name__ == '__main__':
    unittest.main()

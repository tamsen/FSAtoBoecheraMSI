import unittest

from signal_processing import peak_analysis


class MyTestCase(unittest.TestCase):
    def test_something(self):


        # [[91.61385072297661, 1342.1999999999998], [93.81033643218078, 20069.1],
        #        [97.86754348786941, 2004.3000000000002], [99.90795042903845, 18197.400000000005],
        #        [101.91874606301356, 15506.700000000003]]

        #[[93.8, 20069.1], [99.9, 18197.400000000005], [101.9, 15506.700000000003]]
        #[[89.9, 823.1000000000001], [93.9, 13956.1], [99.9, 12466.0], [101.8, 11292.900000000001]]
        typical_stutter = 3.5

        #
        peaks = [[89.9, 823.1000000000001], [92.10233688946954, 855.5999999999999],\
               [93.9, 13956.1], [97.86468582458026, 1325.6], [99.90770391497513, 12466.0],\
               [101.83534326055506, 11292.900000000001], [103.01833113835109, 2365.2]]

        #peaks = [[89.9, 823.1000000000001], [93.9, 13956.1], [99.9, 12466.0], [101.8, 11292.900000000001]]
        step_width_left = 10
        step_proportion_left = 0.4
        step_width_right = 3
        step_proportion_right = 0.4

        while True:
            cleaned_peaks = peak_analysis.step_check(peaks, step_width_left, step_proportion_left,
                                       step_width_right, step_proportion_right)
            if peaks == cleaned_peaks:
                break

        #peaks = cleaned_peaks
        #peaks = peak_analysis.step_check(original_peaks, step_width_left, step_proportion_left,
        #                   step_width_right, step_proportion_right)

        self.assertEqual(peaks[0][0], 93.9)


if __name__ == '__main__':
    unittest.main()

import unittest
import os
import test_globals
import fsa_file_processor
import file_io.xml_file_readers
from signal_processing import peak_analysis


class TestFinalCalls(unittest.TestCase):

    def test_stutter_check_onBF20(self):

        BF20_raw_peaks= [[203.4, 2740.6], [204.3, 1857.1000000000001], [206.3, 559.6999999999999], [208.4, 3363.2000000000007],
     [210.5, 853.4000000000001], [212.5, 2612.7], [213.6, 3730.4], [214.6, 1974.6000000000001], [215.6, 1160.2],
     [216.5, 960.4], [218.5, 2679.9]]

        merge_peaks_closer_than_this = 0.5
        results = peak_analysis.stutter_check_2(BF20_raw_peaks, merge_peaks_closer_than_this,
                                                peak_analysis.take_run_maximum)

        self.assertEqual(results, [[203.4, 2740.6],[204.3, 1857.1000000000001],[206.3, 559.6999999999999],
                        [208.4, 3363.2000000000007], [210.5, 853.4000000000001], [212.5, 2612.7],
                       [213.6, 3730.4],[214.6, 1974.6000000000001],[215.6, 1160.2],
                        [216.5, 960.4], [218.5, 2679.9]])

    def test_stutter_check_2(self):
        merge_peaks_closer_than_this = 2.0
        peaks = [[1.2, 2], [3.1, 7], [4.8, 1], [7, 42], [9.2, 2], [10.1, 7], [14.5, 101], [14.8, 42], [23, 1]]
        results = peak_analysis.stutter_check_2(peaks, merge_peaks_closer_than_this,
                                                peak_analysis.take_run_maximum)
        self.assertEqual(results, [[3.1, 7], [7, 42], [10.1, 7], [14.5, 101], [23, 1]])

        merge_peaks_closer_than_this = 2.0
        peaks = [[1.2, 2], [3.1, 7], [4.8, 1], [7, 42], [9.2, 2], [10.1, 7], [14.5, 101], [14.8, 42], [23, 1]]
        results = peak_analysis.stutter_check_2(peaks, merge_peaks_closer_than_this,
                                                peak_analysis.take_run_centroid)

        self.assertEqual(results,
                         [[3.1, 7], [7.0, 42], [10.1, 7], [14.5, 101],
                          [23.0, 1]])

    def test_PS1_ICE3_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS1_A04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        ice3_results = FSA_file_results.MSI_loci_results_by_loci["ICE3"]
        ice3_raw_alleles_called = ice3_results.raw_alleles_called
        ice3_filtered_alleles_called = ice3_results.filtered_alleles_called
        ice3_final_alleles_called = ice3_results.final_alleles_called

        print("ice3 raw alleles called:" + str(ice3_raw_alleles_called))
        print("filtered alleles called:" + str(ice3_filtered_alleles_called))
        print("final alleles called:" + str(ice3_final_alleles_called))
        # expected are 69,89,131

        self.assertEqual(7, len(ice3_raw_alleles_called))
        self.assertEqual(3, len(ice3_filtered_alleles_called))
        self.assertEqual(3,  len(ice3_final_alleles_called))
        self.assertEqual([69.0, 89.0, 131.0], ice3_final_alleles_called)


    #super-hard one
    def test_PS1_BF20_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS1_A04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["BF20"]
        raw_alleles_called = results.raw_alleles_called
        filtered_alleles_called = results.filtered_alleles_called
        final_alleles_called = results.final_alleles_called

        print("BF20 raw alleles called:" + str(raw_alleles_called))
        print("filtered alleles called:" + str(filtered_alleles_called))
        print("final alleles called:" + str(final_alleles_called))

        self.assertEqual(7, len(raw_alleles_called)) #or 11, if you use smaller Kernel
        self.assertEqual(3, len(filtered_alleles_called))
        self.assertEqual(3, len(final_alleles_called))
        self.assertEqual([209.0,214.0,219.0], final_alleles_called) #close enough. would be better if 214 was 218



    # this is an easy one!!
    def test_PS1_A1_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS1_A04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["A1"]
        raw_alleles_called = results.raw_alleles_called
        filtered_alleles_called = results.filtered_alleles_called
        final_alleles_called = results.final_alleles_called

        print("BF20 raw alleles called:" + str(raw_alleles_called))
        print("filtered alleles called:" + str(filtered_alleles_called))
        print("final alleles called:" + str(final_alleles_called))
        # expected is 238

        self.assertEqual(1, len(raw_alleles_called))
        self.assertEqual(1, len(filtered_alleles_called))
        self.assertEqual(1, len(final_alleles_called))
        self.assertEqual([238.0], final_alleles_called)


    def test_PS2_BF11_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS2_B04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["BF11"]
        raw_alleles_called = results.raw_alleles_called
        filtered_alleles_called = results.filtered_alleles_called
        final_alleles_called = results.final_alleles_called


        print("raw alleles called:" + str(raw_alleles_called))
        print("filtered alleles called:" + str(filtered_alleles_called))
        print("final alleles called:" + str(final_alleles_called))

        self.assertEqual(3, len(raw_alleles_called))
        self.assertEqual(2, len(filtered_alleles_called))
        self.assertEqual(2, len(filtered_alleles_called))
        #self.assertEqual([80,86,95], final_alleles_called)
        self.assertEqual([80, 86], final_alleles_called) #deliberately will not call 95 b/c michael w did not

    def test_PS2_ICE14_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS2_B04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["ICE14"]
        raw_alleles_called = results.raw_alleles_called
        filtered_alleles_called = results.filtered_alleles_called
        final_alleles_called = results.final_alleles_called


        print("raw alleles called:" + str(raw_alleles_called))
        print("filtered alleles called:" + str(filtered_alleles_called))
        print("final alleles called:" + str(final_alleles_called))

        self.assertEqual(3, len(raw_alleles_called))
        self.assertEqual(2, len(filtered_alleles_called))
        self.assertEqual(2, len(filtered_alleles_called))
        self.assertEqual([214,217], final_alleles_called)

    def test_PS2_C8_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS2_B04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["C8"]
        raw_alleles_called = results.raw_alleles_called
        filtered_alleles_called = results.filtered_alleles_called
        final_alleles_called = results.final_alleles_called


        print("raw alleles called:" + str(raw_alleles_called))
        print("filtered alleles called:" + str(filtered_alleles_called))
        print("final alleles called:" + str(final_alleles_called))

        self.assertEqual(3, len(raw_alleles_called))
        self.assertEqual(3, len(filtered_alleles_called))
        self.assertEqual(3, len(filtered_alleles_called))
        self.assertEqual([232,242,253], final_alleles_called) #ideally, wuold be 254 instead of 253



    def test_PS3_BF9_and_BF18_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS3_C04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["BF9"]
        self.assertEqual(8, len(results.raw_alleles_called))
        self.assertEqual(3, len(results.filtered_alleles_called))
        self.assertEqual(3, len(results.final_alleles_called))
        self.assertEqual([82,87,94], results.final_alleles_called) #ideally 82,92, 94

        results = FSA_file_results.MSI_loci_results_by_loci["BF18"]
        self.assertEqual(1, len(results.raw_alleles_called))
        self.assertEqual(1, len(results.filtered_alleles_called))
        self.assertEqual(1, len(results.final_alleles_called))
        self.assertEqual([117], results.final_alleles_called)

    def test_PS3_E9_final_calls(self):

        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                    "BD1200PSE_F04.fsa")
        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["E9"]
        self.assertEqual(5, len(results.raw_alleles_called))
        self.assertEqual(2, len(results.filtered_alleles_called))
        self.assertEqual(2, len(results.final_alleles_called))
        self.assertEqual([197,210], results.final_alleles_called)


    def test_PS4_BF3_and_BF19_and_B6_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS4_D04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["BF3"]
        self.assertEqual(3, len(results.raw_alleles_called))
        self.assertEqual(3, len(results.filtered_alleles_called))
        self.assertEqual(3, len(results.final_alleles_called))
        self.assertEqual([93,99,101], results.final_alleles_called)


        results = FSA_file_results.MSI_loci_results_by_loci["BF19"]
        self.assertEqual(2, len(results.raw_alleles_called))
        self.assertEqual(2, len(results.filtered_alleles_called))
        self.assertEqual(2, len(results.final_alleles_called))
        self.assertEqual([142,150], results.final_alleles_called)


        results = FSA_file_results.MSI_loci_results_by_loci["B6"]
        self.assertEqual(3, len(results.raw_alleles_called))
        self.assertEqual(2, len(results.filtered_alleles_called))
        self.assertEqual(2, len(results.final_alleles_called))

        # deliberately will not call 309/310 b/c michael w did not
        # 309 would be better than 310, but its fine
        self.assertEqual([301,303], results.final_alleles_called)

    def test_PS5_BF15_and_Bdru266_and_A3_final_calls(self):
        panel_info = file_io.xml_file_readers.readPanelXml(test_globals.GLOBAL_panel_file)
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_final_calls",
                                "BD1200PS5_E04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        results = FSA_file_results.MSI_loci_results_by_loci["BF15"]
        self.assertEqual(5, len(results.raw_alleles_called))
        self.assertEqual(3, len(results.filtered_alleles_called))
        self.assertEqual(3, len(results.final_alleles_called))
        self.assertEqual([91,93,99], results.final_alleles_called)


        results = FSA_file_results.MSI_loci_results_by_loci["Bdru266"]
        self.assertEqual(4, len(results.raw_alleles_called))
        self.assertEqual(1, len(results.filtered_alleles_called))
        self.assertEqual(1, len(results.final_alleles_called))
        self.assertEqual([129], results.final_alleles_called)


        results = FSA_file_results.MSI_loci_results_by_loci["A3"]
        self.assertEqual(1, len(results.raw_alleles_called))
        self.assertEqual(1, len(results.filtered_alleles_called))
        self.assertEqual(1, len(results.final_alleles_called))
        self.assertEqual([248], results.final_alleles_called)


if __name__ == '__main__':
    unittest.main()

import unittest
import os
import test_globals
import fsa_file_processor
import file_io.xml_file_readers

class TestAlleleCalls(unittest.TestCase):

    def test_PS1(self):

        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file= os.path.join(test_globals.GLOBAL_test_input_dir,"test_num_allele_calls",
                               "BD1200CNTRL-PS1-A4_H02.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        ice3_results= FSA_file_results.MSI_loci_results_by_loci["ICE3"].raw_alleles_called
        bf20_results= FSA_file_results.MSI_loci_results_by_loci["BF20"].raw_alleles_called
        a1_results= FSA_file_results.MSI_loci_results_by_loci["A1"].raw_alleles_called
        self.assertEqual(8,len(ice3_results))
        self.assertEqual(9,len(bf20_results))
        self.assertEqual(1,len(a1_results))

    def test_PS2(self):

        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file= os.path.join(test_globals.GLOBAL_test_input_dir,"test_num_allele_calls",
                               "BD1200CNTRL-PS2-B4_A03.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        bf11_results = FSA_file_results.MSI_loci_results_by_loci["BF11"].raw_alleles_called
        ice14_results = FSA_file_results.MSI_loci_results_by_loci["ICE14"].raw_alleles_called
        c8_results = FSA_file_results.MSI_loci_results_by_loci["C8"].raw_alleles_called

        self.assertEqual(len(bf11_results), 3)
        self.assertEqual(len(ice14_results), 4)
        self.assertEqual(len(c8_results), 4)

    def test_PS3(self):

        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file= os.path.join(test_globals.GLOBAL_test_input_dir,"test_num_allele_calls",
                               "BD1200CNTRL-PS3-C4_B03.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        bf9_results = FSA_file_results.MSI_loci_results_by_loci["BF9"].raw_alleles_called
        bf18_results = FSA_file_results.MSI_loci_results_by_loci["BF18"].raw_alleles_called
        e9_results = FSA_file_results.MSI_loci_results_by_loci["E9"].raw_alleles_called

        self.assertEqual(len(bf9_results), 8)
        self.assertEqual(len(bf18_results), 2)
        self.assertEqual(len(e9_results), 4)


    def test_PS4(self):

        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file= os.path.join(test_globals.GLOBAL_test_input_dir,"test_num_allele_calls",
                               "BD1200CNTRL-PS4-D4_C03.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        bf3_results = FSA_file_results.MSI_loci_results_by_loci["BF3"].raw_alleles_called
        bf19_results = FSA_file_results.MSI_loci_results_by_loci["BF19"].raw_alleles_called
        b6_results = FSA_file_results.MSI_loci_results_by_loci["B6"].raw_alleles_called

        self.assertEqual(3, len(bf3_results))
        self.assertEqual(3, len(bf19_results))
        self.assertEqual(3, len(b6_results))

    def test_PS5(self):
        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_num_allele_calls",
                                "BD1200PS5c5_G02.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        bf15_results = FSA_file_results.MSI_loci_results_by_loci["BF15"].raw_alleles_called
        bdru266_results = FSA_file_results.MSI_loci_results_by_loci["Bdru266"].raw_alleles_called
        a3_results = FSA_file_results.MSI_loci_results_by_loci["A3"].raw_alleles_called

        self.assertEqual(5, len(bf15_results))
        self.assertEqual(2, len(bdru266_results))
        self.assertEqual(1, len(a3_results))


    def test_PS4_lost_peak(self):

        panel_info = file_io.xml_file_readers.readPanelXml("../data/Panel.xml")
        fsa_file= os.path.join(test_globals.GLOBAL_test_input_dir,"test_num_allele_calls",
                               "BD1200PS4_D04.fsa")

        FSA_file_results = fsa_file_processor.process_fsa_file(
            fsa_file, panel_info, test_globals.GLOBAL_test_output_dir)

        bf3_results = FSA_file_results.MSI_loci_results_by_loci["BF3"].raw_alleles_called
        bf19_results = FSA_file_results.MSI_loci_results_by_loci["BF19"].raw_alleles_called
        b6_results = FSA_file_results.MSI_loci_results_by_loci["B6"].raw_alleles_called

        self.assertEqual(3, len(bf3_results))
        self.assertEqual(2, len(bf19_results))
        self.assertEqual(3, len(b6_results))

if __name__ == '__main__':
    unittest.main()

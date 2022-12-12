import os
import lxml
import requests
import pandas as pd
import test_globals
import unittest

from file_io import query_file


class TestQuery(unittest.TestCase):

    # 'http://sites.biology.duke.edu/windhamlab/files/TD21RP25_SearchResults_ASscore.xls'" >
    # URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Daddy_finder_batch.py"
    # URL3 = "https://sites.biology.duke.edu/windhamlab/files/BD1200CNTRL_SearchResults_ASscore.xls"
    # URL2 = "http://sites.biology.duke.edu/windhamlab/files/TD21RP25_SearchResults_ASscore.xls"

    def setUp(self):

        if not (os.path.exists(test_globals.GLOBAL_test_output_dir)):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

        out_folder = os.path.join(test_globals.GLOBAL_test_output_dir, "test_BMW_query")
        if not (os.path.exists(out_folder)):
            os.makedirs(out_folder)

    def test_html_parsing(self):
        in_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_BMW_query", "BatchQuery_Output.html")
        results = query_file.get_data_from_Search_database_batch_response(in_file)
        print(str(results))

    def test_response(self):

        URL2 = "http://sites.biology.duke.edu/windhamlab/files/TD21RP25_SearchResults_ASscore.xls"
        results2 = query_file.submit_plain_query(URL2)

        out_file = os.path.join(test_globals.GLOBAL_test_output_dir, "test_BMW_query", "TD21RP25_BatchQuery_Ouput.txt")

        with open(out_file, 'wb') as f:
            f.write(results2)

        with open(out_file, 'r') as f:
            lines = f.readlines()

        print("third line:" + lines[2])

    #def test_adjusted_calls(self):
    def adjusted_calls(self): #not a real unit test. Just kept handy here: so I can re-run batchqueries thrugh the BMW

        #in_dir = "/home/tamsen/Data/Eton/by_site/Bennetville/" + \
        #          "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_11_55_23"
        #in_file = os.path.join(in_dir, "BatchQueryForRecall_02_12_2022_11_56_30.csv")

        in_dir = "/home/tamsen/Data/Eton/by_site/DanaGibbs2021/" + \
                  "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_13_03_41"
        in_file = os.path.join(in_dir, "BatchQuery_ForRecall_02_12_2022_13_07_32.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/EmmaLake/" + \
                  "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_21_07"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_24_12.tsv")

        #in_dir = "/home/tamsen/Data/Eton/by_site/RedPeak2021/" + \
        # FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_14_53_07"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_14_57_49.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/GardiskyLake2021/" + \
                    "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_30_36"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_33_03.tsv")


        in_dir = "/home/tamsen/Data/Eton/by_site/GardiskyLake2022/" + \
                    "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_11_56_43"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_11_57_32.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/MonoPass/" + \
                    "FSA_to_microsat_script_results__v1.0.0.6__08_12_2022_12_30_42"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_08_12_2022_12_31_12.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/MountConness/" + \
                    "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_24_51"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_25_01.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/MountLassen/" + \
                    "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_33_25"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_34_15.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/RubyLake2020and2021/" + \
                 "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_11_36"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_19_11.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/SaddlebagLake2021/" + \
                 "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_25_03"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_29_27.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/SaddlebagLake-LakeHelen/" + \
            "FSA_to_microsat_script_results__v1.0.0.6__09_12_2022_09_58_53"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_09_12_2022_10_01_02.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/TwinLakes/" + \
            "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_03_35"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_04_11.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/WalkerLake/" + \
            "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_20_37"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_21_03.tsv")

        in_dir = "/home/tamsen/Data/Eton/by_site/WhiteMountain/" + \
            "FSA_to_microsat_script_results__v1.0.0.6__02_12_2022_12_30_16"
        in_file = os.path.join(in_dir, "BatchQuery_RecallBySample_02_12_2022_12_30_34.tsv")

        URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Search_database_batch.py"
        results1 = query_file.submit_file_query(in_file, URL1)
        out_file = os.path.join(in_dir, "Recall_Output.html")



        with open(out_file, 'wb') as f:
            f.write(results1)

    def test_query(self):

        # cgi - bin / Search_database_byAllele.py
        # Search_database_batch.py
        URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Search_database_batch.py"
        in_file = os.path.join(test_globals.GLOBAL_test_input_dir, "test_BMW_query", "BatchQuery_b.tsv")
        results1 = query_file.submit_file_query(in_file, URL1)
        tables = pd.read_html(results1)  # Returns list of all tables on page

        out_file = os.path.join(test_globals.GLOBAL_test_output_dir, "test_BMW_query", "BatchQuery_Ouput.html")

        with open(out_file, 'wb') as f:
            f.write(results1)

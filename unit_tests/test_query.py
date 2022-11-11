import sys
#import xlwt
#import openpyxl as xl
import requests
#import html5lib

import unittest


class TestQuery(unittest.TestCase):

    def test_query(self):

        URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Daddy_finder_batch.py"
        file = "/home/tamsen/Data/Eton/Frag_Order_16664/FSA_to_microsat_script_results_11_11_2022_11_05_21/" + \
                "BatchQuery_FinalCallsBySample11_11_2022_11_05_58.tsv"

        results1 = self.submit_query(file,URL1)
        print(results1)

        # to download query
        URL2 = "https://sites.biology.duke.edu/windhamlab/files/TD21RP21_SearchResults_ASscore.xls"
        results2 = self.submit_query(file,URL2)
        print(results2)

        file2 = "/home/tamsen/Data/Eton/Frag_Order_16664/FSA_to_microsat_script_results_11_11_2022_11_05_21/" + \
                "BatchQuery_Ouput.tsv"

        with open(file2 , 'wb') as f:
            f.write(results2)

    def submit_query(self,query_filename, URL):
        files = {'file': open(query_filename, 'rb')}
        print("Posting ", query_filename)
        response = requests.post(URL, files=files)
        return response.content

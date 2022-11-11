import os
import requests
import test_globals

import unittest

from file_io import query_file


class TestQuery(unittest.TestCase):


    def setUp(self):

        if not (os.path.exists(test_globals.GLOBAL_test_output_dir)):
            os.makedirs(test_globals.GLOBAL_test_output_dir)

        out_folder=os.path.join(test_globals.GLOBAL_test_output_dir, "test_BMW_query")
        if not (os.path.exists(out_folder)):
            os.makedirs(out_folder)


    def test_response(self):

        URL2 = "http://sites.biology.duke.edu/windhamlab/files/TD21RP25_SearchResults_ASscore.xls"
        results2 = query_file.submit_plain_query(URL2)

        out_file=os.path.join(test_globals.GLOBAL_test_output_dir, "test_BMW_query","TD21RP25_BatchQuery_Ouput.txt")

        with open(out_file, 'wb') as f:
            f.write(results2)

        with open(out_file, 'r') as f:
           lines = f.readlines()

        print("third line:" + lines[2])

    def test_query(self):

        #https: // stackoverflow.com / questions / 11790535 / extracting - data -
        #from
        #-html - table
        #cgi - bin / Search_database_byAllele.py
        #Search_database_batch.py
        #URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Daddy_finder_batch.py"
        URL1 = "https://sites.biology.duke.edu/windhamlab/cgi-bin/Search_database_batch.py"
        in_file=os.path.join(test_globals.GLOBAL_test_input_dir, "test_BMW_query", "BatchQuery.txt")
        results1 = query_file.submit_file_query(in_file,URL1)

        #'http://sites.biology.duke.edu/windhamlab/files/TD21RP25_SearchResults_ASscore.xls'" >
        # to download query
        URL3 = "https://sites.biology.duke.edu/windhamlab/files/BD1200CNTRL_SearchResults_ASscore.xls"
        URL2 = "http://sites.biology.duke.edu/windhamlab/files/TD21RP25_SearchResults_ASscore.xls"
        results2 = query_file.submit_plain_query(URL2)

        out_file=os.path.join(test_globals.GLOBAL_test_output_dir, "test_BMW_query","TD21RP25_BatchQuery_Ouput_TAKE3.txt")

        with open(out_file, 'wb') as f:
            f.write(results1)

    def submit_query(self,query_filename, URL):
        files = {'file': open(query_filename, 'rb')}
        print("Posting ", query_filename)
        response = requests.post(URL, files=files)
        return response.content

    def submit_plain_query(self,URL2):
        response = requests.post(URL2)
        return response.content
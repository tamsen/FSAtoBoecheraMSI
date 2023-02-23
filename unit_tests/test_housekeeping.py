import os
import shutil
import unittest


class test_housekeeping(unittest.TestCase):

    def test_move_and_rename_TD_files(self):

        data_path_in = "/home/tamsen/Data/Eton/Frag_Order_16830/a_raw"
        data_path_out = "/home/tamsen/Data/Eton/Frag_Order_16830/a_renamed"
        all_files = os.listdir(data_path_in)
        for f in all_files:
                print("\tfile:" + f)
                original_file = os.path.join(data_path_in, f)

                if "a-" in f:
                    new_name = f.replace("a-", "-")


                new_file = os.path.join(data_path_out, new_name)
                print("original file: " + original_file)
                print("destination file: " + new_file)
                shutil.copy(original_file, new_file)
                print("copy complete")

        self.assertEqual(True, True)


    def move_and_rename_MW_files(self):

        data_path="/home/tamsen/Data/MWRawFSAData/BMW_Calibration_Batch2"
        all_folders = os.listdir(data_path)
        for d in all_folders:
            print("directory: " + d)
            all_files = os.listdir(os.path.join(data_path,  d))
            for f in all_files:
                print("\tfile:" + f)
                original_file=os.path.join(data_path, d, f)

                if "_QX copy.fsa.txt" in f:
                    new_name=f.replace("_QX copy.fsa.txt",".fsa")
                elif " copy.fsa.txt" in f:
                    new_name = f.replace(" copy.fsa.txt", ".fsa")
                elif "copy.fsa.txt" in f:
                    new_name = f.replace("copy.fsa.txt", ".fsa")

                new_file=os.path.join(data_path, new_name)
                print("original file: " + original_file)
                print("destination file: " + new_file)
                shutil.copy(original_file,new_file)
                print("copy complete")

        self.assertEqual(True, True)

    def move_and_rename_TM_files(self):

        data_path="/home/tamsen/Data/TMRawFSAData/"
        all_folders = os.listdir(data_path)
        for d in all_folders:
            print("directory: " + d)
            all_files = os.listdir(os.path.join(data_path,  d))
            for f in all_files:
                print("\tfile:" + f)
                original_file=os.path.join(data_path, d, f)

                if "_" in f:
                    new_name=f.replace("_","_PS")

                new_file=os.path.join(data_path, new_name)
                print("original file: " + original_file)
                print("destination file: " + new_file)
                shutil.copy(original_file,new_file)
                print("copy complete")

        self.assertEqual(True, True)

if __name__ == '__main__':
    unittest.main()


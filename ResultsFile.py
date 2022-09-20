import os
from datetime import datetime

def WriteResults(outputDir, data):

        resultsFile = os.path.join(outputDir,"Results.csv")
        now = datetime.now()
        day = now.strftime("%d/%m/%Y")
        time = now.strftime("%H:%M:%S")
        data_string = ",".join(data)
        time_stamp_string = ",".join([day,time])

        with open(resultsFile, 'a') as f:
           f.write(time_stamp_string + "," + data_string + "\n")


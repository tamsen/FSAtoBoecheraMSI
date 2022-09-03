import matplotlib
import matplotlib.pyplot as plt

def plotLadder(trace_data_dictionary):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(trace_data_dictionary['DATA105'])
    plt.show()
    plt.savefig("./tmp/plot1" + ".png")
    plt.close()

def plotRawTraceByColor(trace_data_dictionary):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(trace_data_dictionary['DATA105'])
    plt.show()
    plt.savefig("./tmp/plot2" + ".png")
    plt.close()

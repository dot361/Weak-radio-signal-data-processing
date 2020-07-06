
import sys
import os
import re
import numpy as np
import scipy.constants
from astropy.time import Time
import datetime
import matplotlib.pyplot as plt
import glob
from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes
import seaborn as sns
import datetime

#from scipy import signal
#from statsmodels.nonparametric.smoothers_lowess import lowess
#from scipy.signal import savgol_filter
#from pykalman import KalmanFilter
import matplotlib.ticker as ticker



        
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]



dir = sys.argv[1]
logs = list()

for file in sorted(glob.glob(dir + "/*.log")):
    logs.append(file)

logs.sort(key=natural_keys)
totalMeasureTime = datetime.timedelta()
print(logs)
for item in range(0,len(logs)):
    log = LogReaderFactory.getLogReader(LogTypes.SDR, logs[item], "/home/gj/Desktop/test").getLogs()
    scan_1 = log["1s0"]["date"]
    scan_fin = log[str(list(log.keys())[-1])]["date"]
    
    t1 = datetime.datetime.strptime(scan_1, '%Y-%m-%dT%H:%M:%S')
    t2 = datetime.datetime.strptime(scan_fin, '%Y-%m-%dT%H:%M:%S')
    tdiff = t2-t1
    totalMeasureTime = tdiff+totalMeasureTime
    print(totalMeasureTime)

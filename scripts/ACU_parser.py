import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import os
import sys
import json
from collections import deque
import glob
from multiprocessing import Process, Pipe
import time
from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes
import argparse
import scipy.constants
from astropy.time import Time
import csv 
import subprocess
import pandas
import matplotlib.ticker as ticker






logs = LogReaderFactory.getLogReader(LogTypes.SDR, sys.argv[1], "/home/testuser/Desktop/prettylogs").getLogs()
# dir = sys.argv[2]


# q = list()
# for file in sorted(glob.glob(dir + "/*.raw")):
# 	q.append(file)

files = int(logs["header"]["N_scans,N_cal"][0])

plt.figure();
dateList = list()
logAzEl = list()
logEl = list()
print(files)
for index in range(1,files//8):

	pair = ((str(index) + "s0", str(index) + "s1"), (str(index) + "r1", str(index) + "r0"))
	print(pair)
	print(float(logs[pair[0][0]]["AzEl"][1]))
	logAzEl.append(float(logs[pair[0][0]]["AzEl"][1]))
	logAzEl.append(float(logs[pair[0][1]]["AzEl"][1]))
	logAzEl.append(float(logs[pair[1][0]]["AzEl"][1]))
	logAzEl.append(float(logs[pair[1][1]]["AzEl"][1]))

	logEl.append(float(logs[pair[0][0]]["AzEl"][0]))
	logEl.append(float(logs[pair[0][1]]["AzEl"][0]))
	logEl.append(float(logs[pair[1][0]]["AzEl"][0]))
	logEl.append(float(logs[pair[1][1]]["AzEl"][0]))



	#logAzEl.append((float(logs[pair[0][0]]["AzEl"][1]) + float(logs[pair[0][1]]["AzEl"][1]) + float(logs[pair[1][0]]["AzEl"][1]) + float(logs[pair[1][1]]["AzEl"][1])) / 4)
	dateList.append(logs[pair[0][0]]["date"])
	dateList.append(logs[pair[0][1]]["date"])
	dateList.append(logs[pair[1][0]]["date"])
	dateList.append(logs[pair[1][1]]["date"])
	#logEl.append((float(logs[pair[0][0]]["AzEl"][0]) + float(logs[pair[0][1]]["AzEl"][0]) + float(logs[pair[1][0]]["AzEl"][0]) + float(logs[pair[1][1]]["AzEl"][0])) / 4)

	# dateList.append(logs[pair[1][0]]["date"])
	# dateList.append(logs[pair[0][0]]["date"])
	# dateList.append(logs[pair[1][1]]["date"])

print(logAzEl)
print(len(logAzEl))
print(len(dateList))
# plt.subplot(1,2,1)
#plt.plot(AzEl00,label="0,0")
#plt.plot(AzEl01,label="0,1")
#plt.plot(AzEl10,label="1,0")
#plt.plot(AzEl11,label="1,1")


#df = pandas.read_csv(getArgs('aculog'),usecols=[2,3], header=None, sep=' ', skiprows=2)
#print(df[2][0][2])
#splitted = df[2]
acuLogs = list()
for file in sorted(glob.glob(sys.argv[2]+"/STS*")):
    with open (file, "r") as myfile:
    	acuLogs+=myfile.readlines()[1:]
    	myfile.close()

print(dateList);
az_axSts_plst = list()
el_axSts_plst = list()
az_axSts_psoll = list()
el_axSts_psoll = list() 
plotDates = list()
for i in range(0,len(acuLogs),5):
    for j in range(0,len(dateList)):
        if(acuLogs[i].startswith(dateList[j])):
            
            entries = acuLogs[i].split(sep=",")
            
            plotDates.append(entries[0][:-4])
            az_axSts_plst.append(entries[42])
            el_axSts_plst.append(entries[65])
            az_axSts_psoll.append(entries[40])
            el_axSts_psoll.append(entries[63])
    	

# df[2] = df[2].apply(lambda s: ",".join(t.split(" ")[0] for t in s.split(",")))
# elList = list()
# print(df[2][0][3])
# for index in range(1,10000):
# 	elList.append(df[2][index][10])
az_axSts_plst = np.array(az_axSts_plst, dtype=float)
az_axSts_psoll = np.array(az_axSts_psoll, dtype=float)
el_axSts_plst = np.array(el_axSts_plst, dtype=float)
el_axSts_psoll = np.array(el_axSts_psoll, dtype=float)
formattedDates = list()
formattedAcuDates = list()

for	i in dateList:
    formattedDates.append(i[8:])
for	i in plotDates:
    formattedAcuDates.append(i[8:-4])
    

plt.subplot(2,1,1)
plt.plot(formattedDates,logAzEl,label="from_log")
plt.xticks( formattedDates, rotation='vertical')
plt.legend()

plt.plot(formattedAcuDates, el_axSts_plst, label="el_axSts_plst")

plt.plot(formattedAcuDates, el_axSts_psoll, label="el_axSts_psoll")
plt.xticks( formattedAcuDates, rotation='vertical')
plt.legend()

plt.subplot(2,1,2)

plt.plot(formattedDates, logEl, label="from_log")
plt.plot(formattedAcuDates, az_axSts_psoll, label="az_axSts_psoll")
plt.plot(formattedAcuDates, az_axSts_plst, label="az_axSts_plst")
plt.xticks( formattedAcuDates, rotation='vertical')
plt.legend()


#print(az_axSts_plst)
print(el_axSts_plst)
#print(az_axSts_psoll)
print(el_axSts_psoll)
plt.legend()
#plt.yticks(ticks=np.arange(37.405659, 60.989934+1.0, 1.0))


plt.show()
plt.figure()

deltaPsoll = list()
deltaPlst = list()
for i in range(len(el_axSts_plst)):
	deltaPsoll.append(np.abs(logAzEl[i]-el_axSts_psoll[i]))
	deltaPlst.append(np.abs(logAzEl[i]-el_axSts_plst[i]))
deltaPlstSeconds = [i * 3600 for i in deltaPlst]
deltaPsollSeconds = [i * 3600 for i in deltaPsoll]
plt.subplot(2,1,1)

plt.plot(formattedAcuDates, deltaPlstSeconds, label="el_deltaPlst")
plt.plot(formattedAcuDates, deltaPsollSeconds, label="el_deltaPsoll")
plt.xticks( formattedAcuDates, rotation='vertical')
plt.ylabel('Seconds')
plt.legend()
print("avg_AzEl: ")
print( logAzEl)
print("az_axSts_psoll: ")
print(az_axSts_psoll)
print("az-axSts_plst: ")
print(az_axSts_plst)

print("avg_El")
print(logEl)
print("el_axSts_psoll: ")
print(el_axSts_psoll)
print("el_axSts_plst: ")
print(el_axSts_plst)



deltaAzPsoll = list()
deltaAzPlst = list()
for i in range(len(az_axSts_plst)):
	deltaAzPsoll.append(np.abs(logEl[i]-az_axSts_psoll[i]))
	deltaAzPlst.append(np.abs(logEl[i]-az_axSts_plst[i]))
plt.subplot(2,1,2)

deltaAzPlstSeconds = [i * 3600 for i in deltaAzPlst]
deltaAzPsollSeconds = [i * 3600 for i in deltaAzPsoll]

plt.plot(formattedAcuDates, deltaAzPlstSeconds, label="az_deltaPlst")
plt.plot(formattedAcuDates, deltaAzPsollSeconds, label="az_deltaPsoll")
plt.xticks( formattedAcuDates, rotation='vertical')
plt.ylabel('Seconds')
plt.legend()
plt.show()
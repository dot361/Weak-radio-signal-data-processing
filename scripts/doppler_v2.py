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
import pywt
import seaborn as sns
import json
import matplotlib.ticker as ticker
import execnet
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
def call_python_version(Version, Module, Function, ArgumentList):
    gw      = execnet.makegateway("popen//python=python%s" % Version)
    channel = gw.remote_exec("""
        from %s import %s as the_function
        channel.send(the_function(*channel.receive()))
    """ % (Module, Function))
    channel.send(ArgumentList)
    return channel.receive()

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

#velocitys_1 = dopler((freq_1) * (10 ** 6), VelTotal, line_1_F)
def dopler_shift( velocityReceiver, f0):
    c = scipy.constants.speed_of_light
    #velocitySoure = (-((ObservedFrequency / f0) - 1) * c + (velocityReceiver * 1000)) / 1000
    velocitySoure = (f0-(f0/(1-(velocityReceiver/c))))/1000
    #print(velocitySoure)
    return velocitySoure
def dopler(ObservedFrequency, velocityReceiver, f0):
    c = scipy.constants.speed_of_light
    velocitySoure = (-((ObservedFrequency / f0) - 1) * c ) / 1000
    #velocitySoure = ((ObservedFrequency/f0) *c + velocityReceiver*1000)/1000
    return velocitySoure

dir = sys.argv[1]
q = list()
logs = list()
for file in sorted(glob.glob(dir + "/*.dat")):
    q.append(file)

for file in sorted(glob.glob(dir + "/*.log")):
    logs.append(file)
print("measures - ", len(q))
print("logs - ", logs)
q.sort(key=natural_keys)
logs.sort(key=natural_keys)

freq = np.loadtxt(q[0], usecols=(0,), unpack=True)
size = len(freq)
files = len(q)


colors = sns.color_palette("Blues", files)
colors[-1] = "#940071"
sns.set(palette=colors)

source = []
velocities_horizons = []

if any(re.findall(r'panstarrs|atlas', logs[0], re.IGNORECASE)):
#if(("panstarrs" or "atlas") in logs[0]):
    splitted = logs[0].split("_")[0].split("/")
    name = splitted[-1]
    print(name)
    for item in logs:
        print(item)
        log = LogReaderFactory.getLogReader(LogTypes.SDR, item, "/home/gj/Desktop/prettylogs").getLogs()
        #print(log["1s0"]["date"])
        scan_1 = log["1s0"]["date"]
        scan_2 = log["2s0"]["date"]
        result = call_python_version("2.7", "get_ra_dec_horizons", "get_ra_dec",  
                                [scan_1, scan_2, name]) 
        source.append(result[0] + ", " + result[1])
        print(result[2])
        velocities_horizons.append(result[2])
else:
    for item in logs:
        log = LogReaderFactory.getLogReader(LogTypes.SDR, item, "/home/gj/Desktop/prettylogs").getLogs()
        ra_dec = log["header"]["RA,DEC"]
        ra_parsed = re.sub('[hms]', '', ra_dec[0])
        dec_parsed = re.sub('[+dms]', '', ra_dec[1])
        source.append(ra_parsed + ", " + dec_parsed)
        
print("sources provided - ", len(source))
for p in range(1, len(q)+1):
    print(q[p-1])
    if("raw" not in q[p-1]):    
        velocitys_1_avg = np.zeros(int(size/4))
        velocitys_2_avg = np.zeros(int(size/4))
        y_1_left_avg = np.zeros(int(size/4))
        y_2_left_avg = np.zeros(int(size/4))
        y_1_right_avg = np.zeros(int(size/4))
        y_2_right_avg = np.zeros(int(size/4))
    else:
        velocitys_1_avg = np.zeros(int(size/2))
        velocitys_2_avg = np.zeros(int(size/2))
        y_1_left_avg = np.zeros(int(size/2))
        y_2_left_avg = np.zeros(int(size/2))
        y_1_right_avg = np.zeros(int(size/2))
        y_2_right_avg = np.zeros(int(size/2))

    print("=================================================")
    for measure in range(0,p):
        log = LogReaderFactory.getLogReader(LogTypes.SDR, logs[measure], "/home/gj/Desktop/prettylogs").getLogs()
    
        line_2 = "1.6673590"
        line_1 = "1.665402"
        #line = "6.6685192" #6.6685192
        #line = "1.666400000"
        line_1_F = float(line_1) * (10**9)
        line_2_F = float(line_2) * (10**9)
        

        #print("specie", specie, "\n")
        freq = np.loadtxt(q[measure], usecols=(0,), unpack=True)
        left = np.loadtxt(q[measure], usecols=(1,), unpack=True)
        right = np.loadtxt(q[measure], usecols=(2,), unpack=True)

        if("raw" not in q[0]):
            df_div = float(log["header"]["df_div,df"][0])
            #df_div = float(logs["header"]["frst,f0,LO,IF,df_div"][4])

            BW = float(log["header"]["Fs,Ns,RBW"][0])
            f_shift = BW / df_div
            l_spec = len(freq)
            f_step = (freq[l_spec - 1] - freq[0]) / (l_spec - 1)
            n_shift = int(np.rint(f_shift / f_step))
            avg_interval = 0.5  # inner 50%
            si = int(l_spec / 2 - l_spec * avg_interval / 2)
            ei = int(l_spec / 2 + l_spec * avg_interval / 2)
            freq = freq[si:ei]
            left = left[si:ei]
            right = right[si:ei]

        #split the rest freq areas
        freq_1 = freq[:len(freq)//2]    
        freq_2 = freq[len(freq)//2:]    
        left_1 = left[:len(left)//2]    
        left_2 = left[len(left)//2:]
        right_1 = right[:len(right)//2]    
        right_2 = right[len(right)//2:]      

        
        df_1 = dopler_shift( float(velocities_horizons[measure]), line_1_F)
        df_1 = float(log["header"]["f_obs,LO,IF"][0]) - float(log["header"]["v_obs,v_rad,f_lab"][2])
        df_2 = dopler_shift( float(velocities_horizons[measure]), line_2_F)
        df_1 = float(log["header"]["f_obs,LO,IF"][0]) - float(log["header"]["v_obs,v_rad,f_lab"][2])

        print("df1: ", df_1)
        print("df2: ", df_2)        
        #freq_1 = np.asarray([x+df_1 for x in freq_1])
        #freq_2 = np.asarray([x+df_2 for x in freq_2])
        #print(freq_1)
        #print(freq_2)

        #velocitys_1 = dopler((freq_1) * (10 ** 6), float(velocities_horizons[measure]), line_1_F)
        #velocitys_2 = dopler((freq_2) * (10 ** 6), float(velocities_horizons[measure]), line_2_F)
        velocitys_1 = dopler((freq_1) * (10 ** 6), float(log["header"]["v_obs,v_rad,f_lab"][1]), line_1_F)
        velocitys_2 = dopler((freq_2) * (10 ** 6), float(log["header"]["v_obs,v_rad,f_lab"][1]), line_2_F)
        if(p == len(q)):
            print("Plotting plots with freq")
            plt.figure("Individual measurements")
            plt.subplot(2,2,1)
            plt.title("LCP")
            plt.yticks(np.arange(-10, 10, step=0.5))
            plt.xlabel("Frequency (MHz)")
            plt.ylabel("Flux density (Jy)")
            plt.plot(freq_1, left_1)

            plt.grid(True)

            plt.subplot(2,2,2)
            plt.title("RCP")
            plt.xlabel("Frequency (MHz)")
            plt.yticks(np.arange(-10, 10, step=0.5))
            plt.plot(freq_1, right_1)
            plt.grid(True)

            plt.subplot(2,2,3)
            plt.xlabel("Frequency (MHz)")
            plt.ylabel("Flux density (Jy)")
            plt.yticks(np.arange(-10, 10, step=0.5))
            plt.plot(freq_2, left_2)

            plt.grid(True)

            plt.subplot(2,2,4)
            plt.yticks(np.arange(-10, 10, step=0.5))
            plt.xlabel("Frequency (MHz)")
            plt.plot(freq_2, right_2)
            plt.grid(True)




        velocitys_1_avg += velocitys_1
        velocitys_2_avg += velocitys_2
        
        y_1_left_avg += left_1
        y_1_right_avg += right_1
        y_2_left_avg += left_2
        y_2_right_avg += right_2




    print("measurements:", p) 
    #plt.show()
    velocitys_1_avg = velocitys_1_avg/p
    velocitys_2_avg = velocitys_2_avg/p
    
    y_1_left_avg = y_1_left_avg/p
    y_1_right_avg = y_1_right_avg/p
    y_2_left_avg = y_2_left_avg/p
    y_2_right_avg = y_2_right_avg/p
    
    
    plt.figure("1.665402")
    plt.suptitle("1.665402")
    
    plt.subplot(1,2,1)
    plt.yticks(np.arange(-10, 10, step=0.1))
    plt.xticks(np.arange(-300, 300, step=20))
    plt.xlabel("Velocity (km sec$^{-1}$)")
    plt.ylabel("Flux density (Jy)")
    plt.plot(velocitys_1_avg, y_1_left_avg)

    plt.grid(True)
    
    plt.subplot(1,2,2)
    plt.yticks(np.arange(-10, 10, step=0.1))
    plt.xticks(np.arange(-300, 300, step=20))

    plt.xlabel("Velocity (km sec$^{-1}$)")
    plt.plot(velocitys_1_avg, y_1_right_avg)
    plt.grid(True)
    plt.figure("1.6673590")
    plt.suptitle("1.6673590")
    plt.subplot(1,2,1)
    plt.yticks(np.arange(-10, 10, step=0.1))
    plt.xticks(np.arange(-300, 300, step=20))
    plt.xlabel("Velocity (km sec$^{-1}$)")
    plt.ylabel("Flux density (Jy)")
    plt.plot(velocitys_2_avg, y_2_left_avg)
    plt.grid(True)


    plt.subplot(1,2,2)
    plt.yticks(np.arange(-10, 10, step=0.1))
    plt.xticks(np.arange(-300, 300, step=20))
    plt.xlabel("Velocity (km sec$^{-1}$)")
    plt.plot(velocitys_2_avg, y_2_right_avg)
    plt.grid(True)






plt.show()

sns.set()


plt.figure("final 1665")
plt.suptitle("1.665402")
plt.subplot(1,2,1)
plt.yticks(np.arange(-10, 10, step=0.1))
plt.xticks(np.arange(-300, 300, step=20))
plt.xlabel("Velocity (km sec$^{-1}$)")
plt.ylabel("Flux density (Jy)")
plt.plot(velocitys_1_avg, y_1_left_avg)
plt.grid(True)

plt.subplot(1,2,2)
plt.yticks(np.arange(-10, 10, step=0.1))
plt.xticks(np.arange(-300, 300, step=20))

plt.xlabel("Velocity (km sec$^{-1}$)")
plt.plot(velocitys_1_avg, y_1_right_avg)
plt.grid(True)
plt.figure("final 1667")
plt.suptitle("1.6673590")
plt.subplot(1,2,1)
plt.yticks(np.arange(-10, 10, step=0.1))
plt.xticks(np.arange(-300, 300, step=20))
plt.xlabel("Velocity (km sec$^{-1}$)")
plt.ylabel("Flux density (Jy)")
plt.plot(velocitys_2_avg, y_2_left_avg)
plt.grid(True)
plt.subplot(1,2,2)
plt.yticks(np.arange(-10, 10, step=0.1))
plt.xticks(np.arange(-300, 300, step=20))
plt.xlabel("Velocity (km sec$^{-1}$)")
plt.plot(velocitys_2_avg, y_2_right_avg)
plt.grid(True)
plt.show()




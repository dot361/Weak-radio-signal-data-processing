#!/usr/bin/env python3
import numpy as np
from scipy import signal
from scipy import stats
import os
import sys
import json
from collections import deque
import glob
import time
from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes
import argparse
import scipy.constants
from astropy.time import Time
from mpi4py import MPI
import matplotlib.pyplot as plt
import pywt
import warnings



"""
Function applies FFT and Blackman-Harris functions to get the amplitude values for the particular raw data file.
:param file_name: Full path to raw data file
:return: Array of amplitude values.
"""
def read_raw(file_name):
    data = np.fromfile(file_name, np.dtype([('re', np.int16), ('im', np.int16)]))                   #Reads data from file in int16+int16 format
    iq_data = data['re'] + data['im'] * 1j                                                          #Splits data into real and imaginary parts
    iq_data = iq_data * 1e-3                                                                        #Correcting the measurement

    Ns_tot = len(iq_data)
    #The log file must be read before calling the function, as it stores necessary FFT length
    Ns = int(logs["header"]["Fs,Ns,RBW"][1])  #Length of FFT                                                
    Nint = int(Ns_tot / Ns)                   #Amount of bins
    window = signal.blackmanharris(int(Ns))   #Blackman-Harris function coefficients
    spec = np.zeros(int(Ns))                  #Stores end result
    ptr = Ns                                  #Stores the pointer for the end of bin to be processed
    binptr = 0                                #Stores which bin is being processed
    for i in range(Nint):                     #Iterates through all data
        spec_i = np.fft.fft(iq_data[binptr*Ns + ptr-Ns:binptr*Ns + ptr]*window);                    #Applies FFT on data region applying blackman-harris function
        spec = spec + (spec_i.real * spec_i.real + spec_i.imag * spec_i.imag)                       #Combines real and imaginary parts
        #Moves the pointer to fulfill the 66.1% overlapping coefficient
        if(i%2==1):
            ptr = Ns + int(Ns * 0.661)                                                           
        if(i%2==0):
            ptr = Ns + int(Ns * 0.339)
        if(i%3==2 and i!=0):
            #print("next bin")
            binptr = binptr + 1              #Next bin
            ptr = Ns
    spec = np.fft.fftshift(spec)             #Shifts frequency spectre to middle
    spec = spec / Nint;                      #Average calculation

    #An array of amplitude values is produced
    return spec


"""
Function performs data calibration on provided files, based on the Unbiased calibration algorithm.
:param p_sig_left:      Local oscilator signal with noise diode off for left polarization
:param p_sig_right:     Local oscilator signal with noise diode off for right polarization
:param p_sig_on_left:   Local oscilator signal with noise diode on for left polarization
:param p_sig_on_right:  Local oscilator signal with noise diode on for right polarization

:param p_ref_left:      Reference oscilator signal with noise diode off for left polarization
:param p_ref_right:     Reference oscilator signal with noise diode off for right polarization
:param p_ref_on_left:   Reference oscilator signal with noise diode on for left polarization
:param p_ref_on_right:  Reference oscilator signal with noise diode on for right polarization

:param frequencyA:      Array of measurement observed frequencies
:param scan:            Integer of which scan is being processed for fetching information from log files
:param logs:            Observations log file
:param Tsys:            Predicted system temperature values

:return:                Array of flux density values.
"""
def frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, frequencyA, scan, logs, Tsys=0):
    df_div = float(logs["header"]["df_div,df"][0])                                                 #Bandwidth value in new log files
    #df_div = float(logs["header"]["frst,f0,LO,IF,df_div"][4])                                     #Bandwidth value in old log files



    #Pointers for data region to be processed.
    BW = float(logs["header"]["Fs,Ns,RBW"][0])
    f_shift = BW / df_div
    l_spec = len(frequencyA)
    f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
    n_shift = int(np.rint(f_shift / f_step))
    avg_interval = 0.5  
    si = int(l_spec / 2 - l_spec * avg_interval / 2)
    ei = int(l_spec / 2 + l_spec * avg_interval / 2)

    #If function is called and Tsys value is provided, it's assumed that an anomaly is detected in the scan.
    if(Tsys == 0 ):
        #Equations 1.2.2, 1.2.3 for each polarization.
        Tsys_off_1_left = float(logs["header"]["Tcal"][0]) * ((p_ref_on_left + p_ref_left) - np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei])) / (2 * np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei]))
        Tsys_off_2_left = float(logs["header"]["Tcal"][1]) * ((p_sig_on_left + p_sig_left) - np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei])) / (2 * np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei]))
        Tsys_off_1_right = float(logs["header"]["Tcal"][0]) * ((p_ref_on_right + p_ref_right) - np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei])) / (2 * np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei]))
        Tsys_off_2_right = float(logs["header"]["Tcal"][1]) * ((p_sig_on_right + p_sig_right) - np.mean(p_sig_on_right[si:ei] - p_sig_right[si:ei])) / (2 * np.mean(p_sig_on_right[si:ei] - p_sig_right[si:ei]))


    else:
        print(scan, " is using a predicted value!")
        #An empty array is created, then filled with predicted system temperature values.
        Tsys_off_1_left = np.empty(len(frequencyA)); Tsys_off_1_left.fill((Tsys[0][scan])/2)
        Tsys_off_2_left = np.empty(len(frequencyA)); Tsys_off_2_left.fill((Tsys[1][scan])/2)
        Tsys_off_1_right = np.empty(len(frequencyA)); Tsys_off_1_right.fill((Tsys[2][scan])/2)
        Tsys_off_2_right = np.empty(len(frequencyA)); Tsys_off_2_right.fill((Tsys[3][scan])/2)


    #Equations 1.2.4, 1.2.5 for each polarization.
    Ta_1_caloff_left = Tsys_off_1_left * (p_sig_left - p_ref_left) / p_ref_left         
    Ta_1_caloff_right = Tsys_off_1_right * (p_sig_right - p_ref_right) / p_ref_right    

    Ta_2_caloff_left = Tsys_off_2_left * (p_ref_left - p_sig_left) / p_sig_left         
    Ta_2_caloff_right = Tsys_off_2_right * (p_ref_right - p_sig_right) / p_sig_right    


    #Equations 1.2.6, 1.2.7 for each polarization.
    Ta_1_calon_left = (Tsys_off_1_left + float(logs["header"]["Tcal"][0])) * (p_sig_on_left - p_ref_on_left) / p_ref_on_left  
    Ta_1_calon_right = (Tsys_off_1_right + float(logs["header"]["Tcal"][1])) * (p_sig_on_right - p_ref_on_right) / p_ref_on_right  

    Ta_2_calon_left = (Tsys_off_2_left + float(logs["header"]["Tcal"][0])) * (p_ref_on_left - p_sig_on_left) / p_sig_on_left  
    Ta_2_calon_right = (Tsys_off_2_right + float(logs["header"]["Tcal"][1])) * (p_ref_on_right - p_sig_on_right) / p_sig_on_right  

    #Equations 1.2.8, 1.2.9 for each polarization
    Ta_sig_left = (Ta_1_caloff_left + Ta_1_calon_left) / 2
    Ta_sig_right = (Ta_1_caloff_right + Ta_1_calon_right) / 2

    Ta_ref_left = (Ta_2_caloff_left + Ta_2_calon_left) / 2
    Ta_ref_right = (Ta_2_caloff_right + Ta_2_calon_right) / 2

    #Frequency shift 
    Ta_sig_left = np.roll(Ta_sig_left, +n_shift)
    Ta_sig_right = np.roll(Ta_sig_right, +n_shift)

    Ta_ref_left = np.roll(Ta_ref_left, -n_shift)
    Ta_ref_right = np.roll(Ta_ref_right, -n_shift)

    #Equation 1.2.10 for each polarization
    Ta_left = (Ta_sig_left + Ta_ref_left) / 2
    Ta_right = (Ta_sig_right + Ta_ref_right) / 2

    #Fetching the coefficients for elevation calculation polynomial
    pair = ((str(scan) + "s1", str(scan) + "r1"), (str(scan) + "r0", str(scan) + "s0"))
    El = (float(logs[pair[0][0]]["AzEl"][1]) + float(logs[pair[0][1]]["AzEl"][1]) + float(logs[pair[1][0]]["AzEl"][1]) + float(logs[pair[1][1]]["AzEl"][1])) / 4
    G_El = logs["header"]["Elev_poly"]
    G_El = [float(gel) for gel in G_El]
    G_ELtmp = [0, 0, 0]
    G_ELtmp[0] = G_El[2]
    G_ELtmp[1] = G_El[1]
    G_ELtmp[2] = G_El[0]
    G_El = G_ELtmp

    #Calculating the flux density
    Sf_left = Ta_left / (( float(logs["header"]["DPFU"][0])) * np.polyval(G_El, El))
    Sf_right = Ta_right / (( float(logs["header"]["DPFU"][1])) * np.polyval(G_El, El))

    #Cuts out the edges of spectrum, since they are not needed for data processing.
    return (Sf_left[si:ei], Sf_right[si:ei], frequencyA[si:ei])


# There is a known issue in OpenMPI implementations, where idle processes sometimes take up 100% of cpu power.
# To prevent this issue, a solution from https://groups.google.com/forum/#!topic/mpi4py/nArVuMXyyZI is used.
def barrier(comm, tag=0, sleep=0.01):
    size = comm.Get_size()
    if size == 1:
        return
    rank = comm.Get_rank()
    mask = 1
    while mask < size:
        dst = (rank + mask) % size
        src = (rank - mask + size) % size
        req = comm.isend(None, dst, tag)
        while not comm.Iprobe(src, tag):
            time.sleep(sleep)
        comm.recv(None, src, tag)
        req.Wait()
        mask <<= 1

"""
Algorithm uses argparse to fetch provided arguments
"""
def parseArguments():
    parser = argparse.ArgumentParser(description='''Uses MPI protocol to process data in provided directory.''')
    parser.add_argument("logFile", help="Experiment log file name", type=str)
    parser.add_argument("dir",help="directory of data",type=str)
    args = parser.parse_args()
    return args

"""
Helper function to for log parsing
"""
def getArgs(key):
    return str(parseArguments().__dict__[key])


"""
Applies wavelet transformation based input.
:param data:        Any kind of integer or float array
:param wavelet:     Name of wavelet in its short form from pywt docs (ex. rbio1.5)
:param threshold:   Transformation threshold, float value.

:return:            Modified data array, based on inputted parameters.
"""
def wavTransform(data, wavelet, threshold):
    w = pywt.Wavelet(wavelet)
    maxlev = pywt.dwt_max_level(len(data), w.dec_len)
    coeffs = pywt.wavedec(data, wavelet, level=maxlev)
    for i in range (1,len(coeffs)):
        coeffs[i] = pywt.threshold(coeffs[i], threshold*max(coeffs[i]))
    resulting = pywt.waverec(coeffs, wavelet)
    return resulting





#This part gets executed by all processes

comm = MPI.COMM_WORLD           #Used by MPI algorithm for communication with other processes, that are launched by the same command
rank = comm.Get_rank()          #Current process rank

logs = LogReaderFactory.getLogReader(LogTypes.SDR, getArgs("logFile"), "./prettylogs").getLogs()     #Reads the provided log file and parses it, creating a "prettylogs" json file.

#Starting from here, each process has a different job.
#First process is responsible for data calibration and data processing after data is read.
if(rank == 0):
    starttime = time.localtime()                        #Stores start time for logging
    dir = getArgs("dir")                                
    #print(dir)

    #Gets all raw data files in directory, sorts them and puts them into a list
    filelist = list()
    for file in sorted(glob.glob(dir + "/*.raw")):
        filelist.append(file)

    #Ensures that results directory is created
    import errno
    if not os.path.exists( dir + '/results'):
        try:
            os.mkdir( dir + '/results')
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


    files = len(filelist)                                      
    index = 0
    #If there are less than 8 raw data files, it is not possible to perform data calibration, so it is assumed that only dat files should be calibrated.
    if(len(filelist) > 8):

        fileNR = 0      
        freqlist = []                                   #Resulting frequency range

        orderedFiles = list()
        scan_count = logs["header"]["N_scans,N_cal"][0]
        start = float(logs["header"]["f_obs,LO,IF"][2]) - float(logs["header"]["Fs,Ns,RBW"][0])/2               #Calculation of starting frequency for new log files
        #start = float(logs["header"]["frst,f0,LO,IF,df_div"][1]) - float(logs["header"]["Fs,Ns,RBW"][0])/2     #Calculation of starting frequency for old log files
        Ns = int(logs["header"]["Fs,Ns,RBW"][1])                        #Amount of elements in the resulting array
        step = float(logs["header"]["Fs,Ns,RBW"][2])                    #Step size between frequency elements

        #Calculation of frequency range
        for f in range(0, Ns):
            freqlist.append(start)
            start += step

        freq = np.array(freqlist)                                       #Converting to numpy friendly form
        
        df_div = float(logs["header"]["df_div,df"][0])                  #Bandwidth for new log files
        #df_div = float(logs["header"]["frst,f0,LO,IF,df_div"][4])      #Bandwidth for old log files

        #Pointers for data region to be processed
        BW = float(logs["header"]["Fs,Ns,RBW"][0])
        f_shift = BW / df_div
        l_spec = len(freq)
        f_step = (freq[l_spec - 1] - freq[0]) / (l_spec - 1)
        n_shift = int(np.rint(f_shift / f_step))
        avg_interval = 0.5  # inner 50%
        si = int(l_spec / 2 - l_spec * avg_interval / 2)
        ei = int(l_spec / 2 + l_spec * avg_interval / 2)
        print("Amount of files: ", files)

        #Orders the files in a better way, changing the format
        while(filelist):
            measurement = [filelist[4], filelist[5], filelist[0], filelist[1], filelist[6], filelist[7], filelist[2], filelist[3]]
            orderedFiles.extend(measurement)
            filelist= filelist[8:]

        print("Files sorted!")

        #Sends file array to management process
        comm.send(orderedFiles, dest=1, tag=0)
        print("Files sent!")
        #barrier(comm)
        
        #Recieves all data from each process and appends it to total array.
        results = list()
        for process in range(2, comm.size):
            rez = comm.recv(source=process,tag=10001)
            results.append(rez)
        #print("Results: ", results)

        #Saves processed data in an easy to read format.
        np.save(logs["header"]["exp_name"]+"_results.npy",np.asarray(results))

        #Lists for each type of signal
        s1 = list()
        s0 = list()
        r0 = list()
        r1 = list()
        
        #Data from each process is stored as an array of directories, that store file names and corresponding data array. 
        #Since for further data processing it is necessary to provide certain sequences of data, it is sorted in different lists, based on the file name.
        for i in range(0, len(results)):
            for j in range(0, len(results[i])):
                if("s0" in results[i][j]["file"]):
                    s0.append(results[i][j])
                if("s1" in results[i][j]["file"]):
                    s1.append(results[i][j])
                if("r0" in results[i][j]["file"]):
                    r0.append(results[i][j])
                if("r1" in results[i][j]["file"]):
                    r1.append(results[i][j])
        #This seems to be the fastest way to sort data.
        s0 = sorted(s0, key=lambda k: k['file']) 
        s1 = sorted(s1, key=lambda k: k['file']) 
        r0 = sorted(r0, key=lambda k: k['file']) 
        r1 = sorted(r1, key=lambda k: k['file'])

        #To display results from data, a new plot is created.
        plt.figure("raw-data-amplitudes")
        plt.suptitle(logs["header"]["exp_name"] + ' raw data plot')

        freq_lo= [x + float(logs["header"]["f_obs,LO,IF"][1]) for x in freq]
        for i in range(0, len(s0)):
            #Matplotlib allows to select which subplot should be used
            if(i%2 == 0):
                plt.subplot(121)                #Left polarization
            else:
                plt.subplot(122)                #Right polarization
            #plots data in corresponding subplot
            plt.plot(freq_lo, s0[i]["data"])
            plt.plot(freq_lo, r0[i]["data"])
            plt.plot(freq_lo, s1[i]["data"])
            plt.plot(freq_lo, r1[i]["data"])

        #Additional plot information
        plt.subplot(121)
        plt.xlabel("Frequency (MHz) with LO")
        plt.ylabel("Amplitude")
        plt.title("LCP")
        plt.grid(True)
        plt.subplot(122)
        plt.xlabel("Frequency (MHz) with LO")
        plt.ylabel("Amplitude")
        plt.title("RCP")
        plt.grid(True)

        #Saves the plot in the measurement results directory
        plt.savefig(dir + 'results/'+ logs["header"]["exp_name"] +'_raw_plots.png')
        print("RAW plot saved")

        #Calculating average system temperature for each spectrum
        Tcal = [float(logs["header"]["Tcal"][0]), float(logs["header"]["Tcal"][1])]
        Tsys_array = [list(), list(), list(), list()]
        for i in range(0, len(s0), 2):
            #measure = [s0[i]["file"], s0[i+1]["file"], r0[i]["file"], r0[i+1]["file"], s1[i]["file"], s1[i+1]["file"], r1[i]["file"], r1[i+1]["file"]]

            #Eq 1.2.2, 1.2.3
            Tsys_off_1_r = Tcal[0]*((r1[i]["data"] + r0[i]["data"]) - np.mean(r1[i]["data"][si:ei] - r0[i]["data"][si:ei]))/(2*np.mean(r1[i]["data"][si:ei] - r0[i]["data"][si:ei]))
            Tsys_off_1_s = Tcal[1]*((s1[i]["data"] + s0[i]["data"]) - np.mean(s1[i]["data"][si:ei] - s0[i]["data"][si:ei]))/(2*np.mean(s1[i]["data"][si:ei] - s0[i]["data"][si:ei]))
            
            Tsys_off_2_r = Tcal[0]*((r1[i+1]["data"] + r0[i+1]["data"]) - np.mean(r1[i+1]["data"][si:ei] - r0[i+1]["data"][si:ei]))/(2*np.mean(r1[i+1]["data"][si:ei] - r0[i+1]["data"][si:ei]))
            Tsys_off_2_s = Tcal[1]*((s1[i+1]["data"] + s0[i+1]["data"]) - np.mean(s1[i+1]["data"][si:ei] - s0[i+1]["data"][si:ei]))/(2*np.mean(s1[i+1]["data"][si:ei] - s0[i+1]["data"][si:ei]))
            
            Tsys_array[0].append(np.mean(Tsys_off_1_r[si:ei]))
            Tsys_array[1].append(np.mean(Tsys_off_1_s[si:ei]))
            Tsys_array[2].append(np.mean(Tsys_off_2_r[si:ei]))
            Tsys_array[3].append(np.mean(Tsys_off_2_s[si:ei]))


        #System temperature plots
        plt.figure("temperature-plots")
        plt.plot(Tsys_array[0])
        plt.plot(Tsys_array[1])
        plt.plot(Tsys_array[2])
        plt.plot(Tsys_array[3])
        np.save(logs["header"]["exp_name"]+"_s0.npy",np.asarray(Tsys_array[1]))
        np.save(logs["header"]["exp_name"]+"_r0.npy",np.asarray(Tsys_array[0]))
        np.save(logs["header"]["exp_name"]+"_s1.npy",np.asarray(Tsys_array[3]))
        np.save(logs["header"]["exp_name"]+"_r1.npy",np.asarray(Tsys_array[2]))

        plt.savefig(dir + 'results/'+ logs["header"]["exp_name"] +'_tsys_plots.png')

        poly_array = [list(), list(), list(), list()]
        poly_array[0] = wavTransform(Tsys_array[0], 'rbio1.5', 3)
        poly_array[1] = wavTransform(Tsys_array[1], 'rbio1.5', 3)
        poly_array[2] = wavTransform(Tsys_array[2], 'rbio1.5', 3)
        poly_array[3] = wavTransform(Tsys_array[3], 'rbio1.5', 3)

        #Gets the mean value of all system temperature arrays
        Tsys_pol_mean = (np.asarray(Tsys_array[0][:]) + np.asarray(Tsys_array[1][:]) + np.asarray(Tsys_array[2][:]) + np.asarray(Tsys_array[3][:]))/4
        #Calculation of z-score
        z = np.abs(stats.zscore(Tsys_pol_mean)) 
        Tsys_outl_thres = 2;
        #Detects where z-score is above threshold.
        outl = np.where(z > Tsys_outl_thres)
        outlier_count = len(outl[0])
        print("Tsys outliers found: %d" % (outlier_count))
        print(outl)

        

        print("Started calibration!")
        Sf_lefts2 = list()                                 #Resulting flux density spectrum array for left polarization
        Sf_rights2 = list()                                #Resulting flux density spectrum array for right polarization   
        dmged_scans = list()                               #Stores the detected scans for logging

        for i in range(0, len(s0), 2):
            #measure = [s0[i]["file"], s0[i+1]["file"], r0[i]["file"], r0[i+1]["file"], s1[i]["file"], s1[i+1]["file"], r1[i]["file"], r1[i+1]["file"]]
            
            #If anomaly is not detected in the scan, system temperatures are passed to the calibration process.
            if(int((i/2)+1) not in outl[0]):
                Sf_left2, Sf_right2, frequencyA2 = frequency_shifting(s0[i]["data"], s0[i+1]["data"], r0[i]["data"], r0[i+1]["data"], s1[i]["data"], s1[i+1]["data"], r1[i]["data"], r1[i+1]["data"], freq, int((i/2)+1), logs)
            else:
                Sf_left2, Sf_right2, frequencyA2 = frequency_shifting(s0[i]["data"], s0[i+1]["data"], r0[i]["data"], r0[i+1]["data"], s1[i]["data"], s1[i+1]["data"], r1[i]["data"], r1[i+1]["data"], freq, int((i/2)+1), logs, poly_array )
                print("Error detected at ", int((i/2)+1))
                #To help in anomaly resolving, timestamp of error is logged.
                err = {
                    "scan": int((i/2)+1),
                    "timestamp": logs[str(int((i/2)+1)+1)+"s0"]["date"]
                }
                dmged_scans.append(err)
            #Each result is individually stored in the array.

            Sf_lefts2.append(Sf_left2)
            Sf_rights2.append(Sf_right2)

        #Saving results for further use.
        left_polar = np.asarray(Sf_lefts2)
        right_polar = np.asarray(Sf_rights2)
        np.save(logs["header"]["exp_name"]+"_left.npy",left_polar)
        np.save(logs["header"]["exp_name"]+"_right.npy",right_polar)


        #Initializing resulting data array
        Sf_lefts_np2 = np.zeros(len(Sf_lefts2[0]))
        Sf_rights_np2 = np.zeros(len(Sf_rights2[0]))


        #Combines the arrays
        for s in range(0, len(Sf_lefts2)):
            Sf_lefts_np2 = Sf_lefts_np2 + np.array(Sf_lefts2[s])
            Sf_rights_np2 = Sf_rights_np2 + np.array(Sf_rights2[s])

        #Getting the average from all data
        Sf_lefts_np2 /= len(Sf_lefts2)
        Sf_rights_np2 /= len(Sf_rights2)


        #There are encoding issues for numpy and json formats, so it's converted to list
        left = list(Sf_lefts_np2)
        right = list(Sf_rights_np2)
        
        #add local oscilator value to resulting frequency	
        frequencyA2= [x + float(logs["header"]["f_obs,LO,IF"][1]) for x in frequencyA2]




        #Dictionary for json
        avgarray = {}
        avgarray["left"] = left
        avgarray["right"] = right
        avgarray["xaxis"]= list(frequencyA2)
        
        #save to json file
        with open(dir + 'results/'+ logs["header"]["exp_name"] +'_raw.json', 'w') as outfile:
            json.dump(avgarray, outfile) 
        #save to dat type file
        with open(dir + 'results/'+ logs["header"]["exp_name"]+"_raw.dat", "w") as data_file:
            for i in range(len(frequencyA2[si:ei])):
                data_file.write("%f %f %f\n" % (frequencyA2[si:ei][i], left[i], right[i]))
    
        #Stops the logging timer
        endtime = time.localtime()

        #Extra information for logs
        warnings = list()
        if(outlier_count != 0):
            warnings.append("Damaged scans detected!")
        if(files != int(logs["header"]["N_scans,N_cal"][0])*8):
            warnings.append("File count does not match predicted!")


        #Logging data json dictionary data
        rez_logs = {}
        rez_logs["exp_name"]= logs["header"]["exp_name"]
        rez_logs["start_time"]= time.strftime("%Y-%m-%dT%H:%M:%S", starttime)
        rez_logs["end_time"] = time.strftime("%Y-%m-%dT%H:%M:%S", endtime)
        rez_logs["cores"] = comm.size
        rez_logs["file_count"] = files
        rez_logs["scan_count"] = logs["header"]["N_scans,N_cal"][0]
        rez_logs["exp_start"] = logs["1s0"]["date"]
        rez_logs["exp_end"] = logs[str(list(logs.keys())[-1])]["date"]
        rez_logs["anomalies"] = outlier_count
        if(outlier_count != 0):
            rez_logs["errors"] = dmged_scans

        if(len(warnings) != 0):
            rez_logs["warnings"] = warnings
        #saving the logs
        with open(dir + 'results/'+ logs["header"]["exp_name"] +'_calibr_log.json', 'w') as outfile:
            json.dump(rez_logs, outfile, indent=4) 
        #Stops the process.
        barrier(comm)
        #barrier(comm)
    else:
        #Tells the other processes that there are not enough raw data to process, to shut them down.
        comm.send("NO_RAW", dest=1, tag=0)
        #Stops the process.
        barrier(comm)

#Following if statement gets executed for only the management process
if (rank == 1):
    files = comm.recv(source=0,tag=0)           #Fetches file list
    status = [1]*comm.size                      #Checks initialized amount of processes and creates an array for them
    status[0] = 0                               #Calibration process
    status[1] = 0                               #Management process
    status[2] = 0                               #Dat file process
    #print(status)
    allBusy = False
    


    #If there are no raw data in directory, there is no reason to communicate with other threads.
    if(files != "NO_RAW"):
        #To assist in process synchronization and to check connections a message is sent to all threads that are not busy.
        for proc in range(2,len(status)):
            comm.send("INIT", dest=proc,tag=9998)

        #Iterates over file list, sending each file out.
        for i in range(0,len(files)):
            fileSent = False                    #Flag to check file status
            fileName=files.pop(0)               #Pop is used to prevent file processing repetitions.
            while(not fileSent):
                available = np.where(status)    #where function returns tuple of which boolean variables are True, which is useful since 1 and 0 values are used.
                #If all processes are busy, wait for response from other threads.
                if(len(available[0]) == 0):
                    allBusy = True
                    while(allBusy):
                        process = comm.recv( tag=9999)
                        #print(process, " is done")
                        status[process] = 1     #Recv function returns the process ID when assigned to a variable, so it's used to determine which process is free
                        allBusy = False         
                        #break
                        #print("allBusy false")  
                else:
                    #File gets sent to the first available process
                    comm.send(fileName, dest=available[0][0])
                    fileSent=True
                    status[available[0][0]] = 0 #To prevent the process getting called again before data processing is finished.
                    #break
                
        print("All files sent!")
    #Since there is no way for other threads to tell if a process is finished or not, a command is sent to other threads which will trigger process finish.
    for process in range(1,len(status)):
        comm.send("DONE", dest=process, tag=10000)
    #Prevents idle 
    barrier(comm)


#Following if statement gets executed if process is not the management or calibration process.
if(rank !=0 and rank != 1):
    data = list()                       #Stores processed raw data results
    #Since the dat data process will join others to process raw data, it makes sense to put it with the other processes.
    #Data processing is very similar to process rank 0, so only things relevant to only dat file processing are explained.
    if(rank == 2):
        dir = getArgs("dir")
        q = list()
        #Selects all dat files in directory
        for file in sorted(glob.glob(dir + "/*.dat")):
            q.append(file)

        files = len(q)
        index = 0
        fileNR = 0



        Sf_lefts2 = []
        Sf_rights2 = []
        
        freq = np.loadtxt(q[0], usecols=(0,), unpack=True)          #Since the frequency range is already written in the dat files
        #freq = np.array(freqlist)

        df_div = float(logs["header"]["df_div,df"][0])
        #df_div = float(logs["header"]["frst,f0,LO,IF,df_div"][4])

        BW = float(logs["header"]["Fs,Ns,RBW"][0])
        f_shift = BW / df_div
        l_spec = len(freq)
        f_step = (freq[l_spec - 1] - freq[0]) / (l_spec - 1)
        n_shift = int(np.rint(f_shift / f_step))
        avg_interval = 0.5  # inner 50%
        si = int(l_spec / 2 - l_spec * avg_interval / 2)
        ei = int(l_spec / 2 + l_spec * avg_interval / 2)

        plt.figure()
        plt.suptitle(logs["header"]["exp_name"] + " dat data plot")
        freq_lo= [x + float(logs["header"]["f_obs,LO,IF"][1]) for x in freq]
        

        while(q):
            
            index +=1

            #Since in dat format, data is described in format: 
            # Frequency value       Left amplitude value        Right amplitude value
            #there is only 4 files total. They are taken out of a queue, to prevent data processing errors.
            file2 = q.pop(0)
            file4 = q.pop(0)
            file1 = q.pop(0)
            file3 = q.pop(0)

            fileNR=fileNR+8

            #print(file1, file2, file3, file4)
            #Each spectrum is loaded
            p_sig_left =  np.loadtxt(file1, usecols=(1,), unpack=True)
            p_ref_left = np.loadtxt(file2, usecols=(1,), unpack=True)
            p_sig_on_left = np.loadtxt(file3, usecols=(1,), unpack=True)
            p_ref_on_left =np.loadtxt(file4, usecols=(1,), unpack=True)


            p_sig_right = np.loadtxt(file1, usecols=(2,), unpack=True)
            p_ref_right = np.loadtxt(file2, usecols=(2,), unpack=True)
            p_sig_on_right = np.loadtxt(file3, usecols=(2,), unpack=True)
            p_ref_on_right = np.loadtxt(file4, usecols=(2,), unpack=True)

            #shifting values 
            p_sig_right = np.fft.fftshift(p_sig_right)
            p_ref_right = np.fft.fftshift(p_ref_right)
            p_sig_on_right = np.fft.fftshift(p_sig_on_right)
            p_ref_on_right = np.fft.fftshift(p_ref_on_right)

            p_sig_left = np.fft.fftshift(p_sig_left)
            p_ref_left = np.fft.fftshift(p_ref_left)
            p_sig_on_left = np.fft.fftshift(p_sig_on_left)
            p_ref_on_left = np.fft.fftshift(p_ref_on_left)

            #Plotting results
            plt.subplot(121)
            plt.plot(freq_lo, p_sig_left)
            plt.plot(freq_lo, p_ref_left)
            plt.plot(freq_lo, p_sig_on_left)
            plt.plot(freq_lo, p_ref_on_left)

            plt.subplot(122)
            plt.plot(freq_lo, p_sig_right)
            plt.plot(freq_lo, p_ref_right)
            plt.plot(freq_lo, p_sig_on_right)
            plt.plot(freq_lo, p_ref_on_right)

            #calibration
            Sf_left2, Sf_right2, frequencyA2 = frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right,p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, freq, index, logs)
            Sf_lefts2.append(Sf_left2)
            Sf_rights2.append(Sf_right2)

        plt.subplot(121)
        plt.title('LCP')
        plt.xlabel("Frequency (MHz) with LO")
        plt.ylabel("Amplitude")
        plt.grid(True)
        plt.subplot(122)
        plt.xlabel("Frequency (MHz) with LO")
        #plt.ylabel("Amplitude")
        plt.title('RCP')
        plt.grid(True)

        plt.savefig(dir + 'results/'+ logs["header"]["exp_name"] +'_dat_plots.png')

        #print(Sf_lefts2)
        left_polar = np.asarray(Sf_lefts2)
        right_polar = np.asarray(Sf_rights2)

        Sf_rights_np2 = np.zeros(len(Sf_rights2[0]))
        Sf_lefts_np2 = np.zeros(len(Sf_lefts2[0]))
        errs = 0
        #Because of system temperature calculation difficulities, currently anomaly detection is not implemented for dat format.
        #Resulting calibration data mean value is checked for failed scan detection, detected scans are not taken into the total result.
        for s in range(0, len(Sf_lefts2)):
            if(np.abs(np.mean(Sf_lefts2[s])) < 2 and np.abs(np.mean(Sf_rights2[s])) < 2):
                Sf_lefts_np2 = Sf_lefts_np2 + np.array(Sf_lefts2[s])
                Sf_rights_np2 = Sf_rights_np2 + np.array(Sf_rights2[s])
            else:
                errs = errs+1

        #average
        Sf_lefts_np2 = Sf_lefts_np2 / (len(Sf_lefts2) -errs )
        Sf_rights_np2 = Sf_rights_np2 / (len(Sf_rights2) -errs )

        left = list(Sf_lefts_np2)
        right = list(Sf_rights_np2)
        
        #add local oscilator value to resulting frequency	
        frequencyA2= [x + float(logs["header"]["f_obs,LO,IF"][1]) for x in frequencyA2]

        #dictionary for json
        plt.figure()
        avgarray = {}
        avgarray["left"] = left
        avgarray["right"] = right
        avgarray["xaxis"]= list(frequencyA2)
        plt.subplot(121)
        plt.plot(frequencyA2,left)
        plt.ylabel("Flux density (Jy)")
        plt.xlabel("Frequency (MHz)")
        plt.grid(True)
        plt.subplot(122)
        plt.xlabel("Frequency (MHz)")
        plt.plot(frequencyA2,right)
        plt.grid(True)
        #plt.show()
        plt.savefig(dir + '/results/'+ logs["header"]["exp_name"] +'_dat.png')


        #save to json file

        print("Saving result to: ", dir + '/results/'+ logs["header"]["exp_name"] +'_dat.json')
        with open(dir + '/results/'+ logs["header"]["exp_name"] +'_dat.json', 'w') as outfile:
            json.dump(avgarray, outfile) 



        with open(dir + '/results/'+ logs["header"]["exp_name"]+"_dat.dat", "w") as data_file:
            for i in range(len(frequencyA2)):
                data_file.write("%f %f %f\n" % (frequencyA2[i], left[i], right[i]))
    


        print("DAT calibration done!")
        comm.send(comm.Get_rank(), dest=1, tag=9999)
        #Starts raw data processing.

        #barrier(comm)
    
    #Each raw data process keeps working until an end signal is recieved. 
    #While using while True loops is usually not the best practice, there are technical issues with other types of loop conditions, 
    # since the end condition is somewhat similar to the promise system in some of the web app languages.
    while(True):
        #Gets the message from management process
        file = comm.recv(source=1)

        #Checks the communication between processes.
        if(file == "INIT"):
            print(comm.Get_rank(), " initalized!")
            comm.send(comm.Get_rank(), dest=1, tag=9999)
        #Way to establish a way for the process to end its work, since there is no way for it to know about the amount of remaining work.
        elif(file == "DONE"):
            comm.send(data, dest=0, tag=10001)
            break
        #An actual file name gets sent
        else:   
            #For control
            print(comm.Get_rank(), "recieved", file)
            
            raw = read_raw(file)
            #saves result in a dictionary
            result = {
                "file": file,
                "data": raw
            }
            #Each process stores all of its processed data.
            data.append(result)
            #When data processing is done, a message is sent to management process.
            comm.send(comm.Get_rank(), dest=1, tag=9999)

    #Process ends its lifespan.        
    print(comm.Get_rank(), " recieved DONE.")
    barrier(comm)




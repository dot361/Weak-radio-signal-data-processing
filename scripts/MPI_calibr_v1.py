#!/usr/bin/env python3
import numpy as np
from scipy import signal
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
import csv 
from mpi4py import MPI


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

:return:                Array of flux density values.
"""
def frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right, p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, frequencyA, scan, logs):
    df_div = float(logs["header"]["df_div,df"][0])                                                 #Bandwidth value in new log files
    #df_div = float(logs["header"]["frst,f0,LO,IF,df_div"][4])                                     #Bandwidth value in old log files

    #Pointers for data region to be processed.
    BW = float(logs["header"]["Fs,Ns,RBW"][0])
    f_shift = BW / df_div
    l_spec = len(frequencyA)
    f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
    n_shift = int(np.rint(f_shift / f_step))
    avg_interval = 0.5  # inner 50%
    si = int(l_spec / 2 - l_spec * avg_interval / 2)
    ei = int(l_spec / 2 + l_spec * avg_interval / 2)

    #Equations 1.2.2, 1.2.3
    Tsys_off_1_left = float(logs["header"]["Tcal"][0]) * ((p_ref_on_left + p_ref_left) - np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei])) / (2 * np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei]))
    Tsys_off_2_left = float(logs["header"]["Tcal"][1]) * ((p_sig_on_left + p_sig_left) - np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei])) / (2 * np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei]))

    Tsys_off_1_right = float(logs["header"]["Tcal"][0]) * ((p_ref_on_right + p_ref_right) - np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei])) / (2 * np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei]))
    Tsys_off_2_right = float(logs["header"]["Tcal"][1]) * ((p_sig_on_right + p_sig_right) - np.mean(p_sig_on_right[si:ei] - p_sig_right[si:ei])) / (2 * np.mean(p_sig_on_right[si:ei] - p_sig_right[si:ei]))

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


#This part gets executed by all processes

comm = MPI.COMM_WORLD           #Used by MPI algorithm for communication with other processes, that are launched by the same command
rank = comm.Get_rank()          #Current process rank

logs = LogReaderFactory.getLogReader(LogTypes.SDR, getArgs("logFile"), "./prettylogs").getLogs()     #Reads the provided log file and parses it, creating a "prettylogs" json file.

#Starting from here, each process has a different job.
#First process is responsible for data calibration and data processing after data is read.
if(rank == 0):
    dir = getArgs("dir")
    q = list()
    
    #Gets all raw data files in directory, sorts them and puts them into a list
    for file in sorted(glob.glob(dir + "/*.raw")):
        q.append(file)
    files = len(q)
    index = 0
    #If there are less than 8 raw data files, it is not possible to perform data calibration.
    if(len(q) > 8):

        fileNR = 0


        Sf_lefts2 = []
        Sf_rights2 = []
        freqlist = []

        #Calculate starting frequency from log file
        start = float(logs["header"]["f_obs,LO,IF"][2]) - float(logs["header"]["Fs,Ns,RBW"][0])/2
        #for older log files		
        #start = float(logs["header"]["frst,f0,LO,IF,df_div"][1]) - float(logs["header"]["Fs,Ns,RBW"][0])/2
        Ns = int(logs["header"]["Fs,Ns,RBW"][1])
        step = float(logs["header"]["Fs,Ns,RBW"][2])

        #Calculating frequency range
        for f in range(0, Ns):
            freqlist.append(start)
            start += step

        #To numpy friendly format
        freq = np.array(freqlist)
        print("Amount of files: ", files)
        #Iterates through queue, popping 8 measurements each time.
        while(q):
			
            index +=1
            

            #Files get sent to predetermined processes
            file3 = q.pop(0)	#r0 polar1
            comm.send(file3, dest=3)
            file4 = q.pop(0)	#r0 polar2
            comm.send(file4, dest=4)
            file7 = q.pop(0)	#r1 polar1
            comm.send(file7, dest=7)
            file8 = q.pop(0)	#r1 polar2
            comm.send(file8, dest=8)
            file1 = q.pop(0)	#s0 polar1
            comm.send(file1, dest=1)
            file2 = q.pop(0)	#s0 polar2
            comm.send(file2, dest=2)
            file5 = q.pop(0)	#s1 polar1
            comm.send(file5, dest=5)
            file6 = q.pop(0)	#s1 polar2
            comm.send(file6, dest=6)


            fileNR=fileNR+8

            #Fetches back the data read by other processes
            p_sig_left = comm.recv(source=1)
            p_ref_left = comm.recv(source=3)
            p_sig_on_left = comm.recv(source=5)
            p_ref_on_left = comm.recv(source=7)



            p_sig_right = comm.recv(source=2)
            p_ref_right = comm.recv(source=4)
            p_sig_on_right = comm.recv(source=6)
            p_ref_on_right = comm.recv(source=8)
            
            
            #Calibration process
            Sf_left2, Sf_right2, frequencyA2 = frequency_shifting(p_sig_left, p_sig_right, p_ref_left, p_ref_right,p_sig_on_left, p_sig_on_right, p_ref_on_left, p_ref_on_right, freq, index, logs)

            #Saves the result
            Sf_lefts2.append(Sf_left2)
            Sf_rights2.append(Sf_right2)

		#Sends a message about ending processes
        comm.send("Finished", dest=1)
        comm.send("Finished", dest=2)
        comm.send("Finished", dest=3)
        comm.send("Finished", dest=4)
        comm.send("Finished", dest=5)
        comm.send("Finished", dest=6)
        comm.send("Finished", dest=7)
        comm.send("Finished", dest=8)

        #Save results for further use
        left_polar = np.asarray(Sf_lefts2)
        right_polar = np.asarray(Sf_rights2)
        np.save(logs["header"]["exp_name"]+"_left.npy",left_polar)
        np.save(logs["header"]["exp_name"]+"_right.npy",right_polar)
        Sf_lefts_np2 = np.zeros(len(Sf_lefts2[0]))
        for s in Sf_lefts2:
            Sf_lefts_np2 = Sf_lefts_np2 + np.array(s)
        Sf_rights_np2 = np.zeros(len(Sf_rights2[0]))
        for s in Sf_rights2:
            Sf_rights_np2 = Sf_rights_np2 + np.array(s)
		

        left = list(Sf_lefts_np2)
        right = list(Sf_rights_np2)
        
        #add local oscilator value to resulting frequency	
        frequencyA2= [x + float(logs["header"]["f_obs,LO,IF"][1]) for x in frequencyA2]

        #dict for json
        avgarray = {}
        avgarray["left"] = left
        avgarray["right"] = right
        avgarray["xaxis"]= list(frequencyA2)
        
        #save to json file
        with open(sys.argv[2] + '/results/calibr_rez.json', 'w') as outfile:
            json.dump(avgarray, outfile) 
		# p.kill() 
	
		
#If process isn't process 0, this if statement gets executed.
if(rank!=0):
	while(True):
        #waits for a message from process 0
		file = comm.recv(source=0)
		print(file)
		if (file == "Finished"):
			break
		read_data = read_raw(file)
        #sends back the result
		comm.send(read_data, dest=0)


import urllib
from datetime import datetime,timedelta
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from itertools import cycle


def get_ra_dec(start_time, stop_time, obj):
    lat = 21.847222;
    lon = 57.5593055;
    alt = 10;

    step_size = '10 m'
    if(obj == "panstarrs"):
        obj_name = 'C/2017 T2'
    if(obj == "atlas"):
        #obj_name = 'C/2019 Y4'
        obj_name = '90004453'
    if(obj == "swan"):
        #obj_name = 'C/2019 Y4'
        obj_name = 'C/2020 F8'


    coord_str = str(lat)+','+str(lon)+','+str(alt)
    url = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='"+ obj_name +"'&CENTER='coord'&SITE_COORD='"+ coord_str +"'&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'&START_TIME='"+ start_time +"'&STOP_TIME='" + stop_time+ "'&STEP_SIZE='"+ step_size +"'&QUANTITIES='1,4'&CSV_FORMAT='YES'"
    print(url)
    s = urllib.urlopen(url).read()
    result = ((s.split("$$SOE"))[1].split("$$EOE")[0]).split('\n')
    #print("result", result)
    
    dates = list()
    el = list()
    for i in range(1, len(result)-1):

        result[i] = result[i].replace(" ", "")
        split_rez = result[i].split(',')
        #print(split_rez)
        #print(split_rez[3], split_rez[6])
        #dates.append(split_rez[0])
        dates.append(split_rez[0][:-5] + " \n " + split_rez[0][-5:])

        #print(split_rez[6], float(split_rez[6]))
        el.append(float(split_rez[6]))
        #split_rez[4] = split_rez[4].replace("+","")
    #print(np)
    el = np.array(el)
#    print(el[:] < 15 )
    #print(dates)
    #print(el.dtype)
    el_masked = np.ma.masked_where(el < 15, el, copy=False)
    print(el_masked)
    plt.ylabel("Elevation")
    plt.grid(True)
    axes = plt.gca()
    ticks_to_use = dates[::32*3]
    plt.yticks(np.arange(-360, 360, step=2))

    plt.title("Elevation of " + 'ATLAS (C/2019 Y4)' + " from " + sys.argv[1] + " to " +sys.argv[2])
    plt.plot(dates, el,'k')
    plt.plot(el_masked, 'r', linewidth=2)
    plt.axhline(15,color='k', linestyle='--')

    axes.set_xticks(ticks_to_use)
    plt.xticks(rotation=70)
    plt.show()
    #print(dates)
    #print(split_rez[3],  split_rez[6])


    return split_rez[3], split_rez[4], split_rez[6]	
	
   
   
get_ra_dec(sys.argv[1], sys.argv[2], sys.argv[3])
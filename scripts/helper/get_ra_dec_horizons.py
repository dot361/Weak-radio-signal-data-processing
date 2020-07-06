import urllib
from datetime import datetime,timedelta
import sys

def get_ra_dec(start_time, stop_time, obj):
    lat = 21.847222;
    lon = 57.5593055;
    alt = 10;

    step_size = '5 m'
    if(obj == "panstarrs"):
        obj_name = 'C/2017 T2'
    if(obj == "atlas"):
        #obj_name = 'C/2019 Y4'
        obj_name = '90004453'
    if(obj == "swan"):
        #obj_name = 'C/2019 Y4'
        obj_name = 'C/2020 F8'


    coord_str = str(lat)+','+str(lon)+','+str(alt)
    url = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='"+ obj_name +"'&CENTER='coord'&SITE_COORD='"+ coord_str +"'&MAKE_EPHEM='YES'&TABLE_TYPE='OBSERVER'&START_TIME='"+ start_time +"'&STOP_TIME='" + stop_time+ "'&STEP_SIZE='"+ step_size +"'&QUANTITIES='1,20'&CSV_FORMAT='YES'"
    print(url)
    s = urllib.urlopen(url).read()
    result = ((s.split("$$SOE"))[1].split("$$EOE")[0]).split('\n')
    print("result", result)

    result[1] = result[1].replace(" ", "")
    split_rez = result[1].split(',')
    split_rez[4] = split_rez[4].replace("+","")

    print(split_rez[0], split_rez[3], split_rez[4], split_rez[6])
    return split_rez[3], split_rez[4], split_rez[6]	
	
   
   
#get_ra_dec(sys.argv[1], sys.argv[2], sys.argv[3])

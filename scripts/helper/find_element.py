
import pandas
import sys

with open (sys.argv[1], "r") as myfile:
    data=myfile.readlines()
entries = data[0].split(sep=",")
print(entries)
print(entries.index("el_axSts.p_Soll"))
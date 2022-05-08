#!/usr/bin/python
# Version 1.7
import sys
import os
import datetime
from datetime import datetime
import time
dt_start = datetime.now()
path = '/work/'
dir_list = os.listdir(path) 
print("Files and directories in '", path, "' :") 
  
# prints all files 
#print(dir_list) 
print("Starting Program to Run Version 5")
print(sys.version)    
os.system("java -jar /work/jython.jar /work/GeneralJMONSELdetector.py")
dt_end = datetime.now()
date_diff = (dt_end - dt_start).total_seconds()
print("Difference in seconds is " + str(date_diff))
print("Difference in minutes is " + str(date_diff/60.00))
print("Start Date:  " + str(dt_start.strftime("%Y-%m-%d %H:%M %Z")))
print("End Date:  " + str(dt_end.strftime("%Y-%m-%d %H:%M %Z")))

print("Finished")
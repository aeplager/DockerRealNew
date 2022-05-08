#!/usr/bin/python

import sys
import os
path = '/work/'
dir_list = os.listdir(path) 
print("Files and directories in '", path, "' :") 
  
# prints all files 
#print(dir_list) 
print("hello there! Papa 5")
print(sys.version)  
  
os.system("java -jar /work/jython.jar /work/GeneralJMONSEL.py")
print("Finished")
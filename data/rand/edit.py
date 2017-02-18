import sys,os
from glob import glob

times = {}
for filename in glob("*st"):
   #print "Open", sys.argv[1]

   a = filename.split("_")
   f = open(filename,"w")
   f.write("")
   f.close()

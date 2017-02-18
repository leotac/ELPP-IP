import sys,os
from glob import glob

times = {}
for filename in glob("*elpp"):
   #print "Open", sys.argv[1]

   blocks=filename.split("_")
   inst = blocks[3]
   s = blocks[6]
   t = blocks[7]

   print inst, s, t
   idx = times.setdefault(inst,1)
   os.system("mv " + filename + " " + inst + "_" + str(idx) + "_" + s + "_" + t)
   times[inst] = times[inst] + 1

import sys,os, re
from glob import glob

times = {}
for filename in glob("*st"):
   #print "Open", sys.argv[1]

   a = filter(None, re.split("[_,.]+", filename))
   s = int(a[2])
   t = int(a[3])

   f = open(filename,"w")
   f.write(str(s) + " " + str(t) + "\n")
   f.close()

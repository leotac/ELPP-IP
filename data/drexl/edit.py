import sys,os
from glob import glob

times = {}
for filename in glob("*st"):
   #print "Open", sys.argv[1]

   a = filename.split("_")
   if a[0] == "rnd-d":
      t = int(a[1])
   else:
      t = int(a[1]) + 2

   f = open(filename,"w")
   f.write("1 " + str(t) + "\n")
   f.close()

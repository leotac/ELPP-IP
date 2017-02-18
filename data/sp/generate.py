import sys,os,random
from glob import glob

times = {}
for filename in glob("*st"):

   N = 30 #number of pairs to generate

   a = filename.split("_")
   n = 3353  #number of nodes
   f = open(filename,"w")
   i = 0
   random.seed(2)
   while i < N:
      s = random.randint(1,n)
      t = random.randint(1,n)
      if s != t:
         i += 1
         f.write(str(s) + " " + str(t) + "\n")
   f.close()

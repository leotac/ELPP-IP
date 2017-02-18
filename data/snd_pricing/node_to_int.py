import sys, os

print "Open", sys.argv[1]
filename=sys.argv[1]
f = open(filename)
idxfile = open(filename+".idx","w")
dictfile = open(filename+".dict","w")

n,m = f.readline().strip().split()
n = int(n)
m = int(m)
print n,m
idxfile.write(str(n)+" "+str(m)+"\n")

s,t = f.readline().strip().split()

names = {}
index = {}
for i in xrange(n):
   node = f.readline().strip()
   names[i]=node
   index[node]=i
   idxfile.write(str(i)+"\n")
   dictfile.write(str(i)+" "+node+"\n")
  
for i in xrange(m):
   u,v,c = f.readline().strip().split()
   idxfile.write(str(index[u])+" "+str(index[v])+" "+c+"\n")

f.close()
idxfile.close()
dictfile.close()

idxfilename = filename+".idx"
newfilename = idxfilename.replace(s, str(index[s]))
newfilename = newfilename.replace(t, str(index[t]))

os.system("mv " + idxfilename + " " + newfilename)

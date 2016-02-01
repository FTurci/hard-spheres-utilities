#!/usr/bin/env python
import os
import scipy as sp
from scipy import spatial,stats
import numpy as np
from numpy import random
import sys
import time

start_time = time.time()
np.seterr(divide='ignore', invalid='ignore')
filename= ''

print("#----------------------------------------------------------------")
print("#-----------Tool for determination of G(r)-----------------------")
print("#----------------------------------------------------------------")

if len(sys.argv)<5:
    print("ERROR: missing mandatory argument!")
    print("usage: Gr coordfile  binsize[pixel] xres[micron] zres[micron]")
    print("example SeriesXX_coords.txt 0.5 0.29 0.25")
    exit(11)
if len(sys.argv) >= 5:
  filename= sys.argv[1]
  dr = float(sys.argv[2])   
  xyres=float(sys.argv[3]) 
  zres=float(sys.argv[4])


print("# file %s"%filename)
C= np.genfromtxt(os.path.join(os.getcwd(),filename),skiprows=2,usecols=(1,2,3))


print("# size before cutting borders: ")

print("# x size: %d %d "%(C[:,0].min(),C[:,0].max()))
print("# y size: %d %d "%(C[:,1].min(),C[:,1].max()))
print("# z size: %d %d "%(C[:,2].min(),C[:,2].max()))

border=10

C = C[np.logical_not(np.logical_or(C[:,0]<C[:,0].min()+border, C[:,0]>C[:,0].max()-border))]
C = C[np.logical_not(np.logical_or(C[:,1]<C[:,1].min()+border, C[:,1]>C[:,1].max()-border))]
C = C[np.logical_not(np.logical_or(C[:,2]<C[:,2].min()+150, C[:,2]>C[:,2].max()-150))]


num_particles=len(C)

print("# size after cutting borders: ")
print("# x size: %d %d "%(C[:,0].min(),C[:,0].max()))
print("# y size: %d %d "%(C[:,1].min(),C[:,1].max()))
print("# z size: %d %d "%(C[:,2].min(),C[:,2].max()))


print("# Number of Particles: %d "%num_particles) 
print("# init done %1.3f s"%(time.time() - start_time))
start_time = time.time()

#create random numbers in the same range as the real particles 
ID0=np.random.randint(C[:,0].min(),C[:,0].max()+1,len(C))
ID1=np.random.randint(C[:,1].min(),C[:,1].max()+1,len(C))
ID2=np.random.randint(C[:,2].min(),C[:,2].max()+1,len(C))


ID=np.vstack((ID0,ID1))
ID=(np.vstack((ID,ID2))).T

#take different zres into account
mul=np.array([1,1,xyres/zres])
ID=np.multiply(ID, mul)
#calculate distances 
RID=sp.spatial.distance.pdist(ID, 'euclidean').flatten()

print("# ideal gas done %1.3f s"%(time.time() - start_time))
start_time = time.time()

mul=np.array([1,1,xyres/zres])
C=np.multiply(C, mul)
RC=sp.spatial.distance.pdist(C, 'euclidean').flatten()

print("# particles done %1.3f s"%(time.time() - start_time))
start_time = time.time()

bins=np.arange(1,200,dr)

print("# bin size: %f array [%f,%f]"%(dr,1,RC.max()))
print("# number of bins: %d"%(len(bins)))

import numpy as np
# H,bins,binnumbers=sp.stats.binned_statistic(RC, RC, statistic='count', bins=bins)
# HID,binsID,binnumbersID=sp.stats.binned_statistic(RID, RID, statistic='count', bins=bins)

H,bins=np.histogram(RC, bins=bins)
HID,binsID=np.histogram(RID,  bins=bins)

#calculate bin centers and divide ParticleHist/IdealGasHist
bincenters = 0.5*(bins[1:]+bins[:-1])
hist=H/HID.astype(float)
#take care of 0/0 and x/0, set to 0
hist[np.isnan(hist)] = 0
hist[np.isinf(hist)] = 0
#save result in file
f=open(filename+"_CC_r.hist",'wb')
np.savetxt(f, np.column_stack((bincenters,hist,H,HID)), fmt='%f')

f.close()

print("# binning done %1.3f s"%(time.time() - start_time))
import pylab as pl
pl.plot(bincenters[:len(bincenters)/3],hist[:len(bincenters)/3],'-')
pl.ylim(0,3)
pl.show()
exit(0)

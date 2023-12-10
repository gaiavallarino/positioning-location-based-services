import pandas as pd
import numpy as np
from pyproj import Proj
import matplotlib.pyplot as plt
import math as mt
# import RTK points from csv file (trajectories)
rtk=pd.read_csv('RTK1.csv',sep=',')
rtkx=rtk.x
rtky=rtk.y
# import smartphone points from csv file
smart1=pd.read_csv('smart1.csv',sep=',')
smart2=pd.read_csv('smart2.csv',sep=',')
#convert from geod to utm
p = Proj(proj='utm',zone=32 ,ellps='WGS84', preserve_units=False)
s1x,s1y=p(smart1['longitude'],smart1['latitude'])
s2x,s2y=p(smart2['longitude'],smart2['latitude'])
# plot 2D trajectories
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(rtk['x'], rtk['y'], '.b', label='rtk trajectory')
ax.plot(s1x, s1y, '.g', label='smart1 trajectory')
ax.plot(s2x, s2y, '.r', label='smart2 trajectory')
plt.title("Trajectories")
leg = ax.legend();
# cycle on smartphones points
dist1=[0]*len(rtk)
dist2=[0]*len(rtk)
err1=[0]*len(smart1)
err2=[0]*len(smart2)
#error evaluation for smartphone 1
for i in range(len(smart1)):
	#find the nearest RTK points
    for j in range(len(rtk)):
         dx1=s1x[i]-rtkx[j]
         dy1=s1y[i]-rtky[j]
         dist1[j]=mt.sqrt(dx1**2+dy1**2)        
    minpos1=dist1.index(min(dist1))
    err1[i]=dist1[minpos1]
#error evaluation for smartphone 2    
for i in range(len(smart2)):
    for j in range(len(rtk)):
        dx2=s2x[i]-rtkx[j]
        dy2=s2y[i]-rtky[j]
        dist2[j]=mt.sqrt(dx2**2+dy2**2)
    minpos2=dist2.index(min(dist2))
    err2[i]=dist2[minpos2]
#plot
fig, ax2 = plt.subplots(figsize=(10,6))
ax2.plot(err1, '.g', label='errors for smartphone 1')
ax2.plot(err2, '.r', label='errors for smartphone 2')
plt.title("errors")
leg = ax2.legend();
#statistics
print('Mean error for smartphone 1: ',np.mean(err1))
print('Mean error for smartphone 2: ',np.mean(err2))
print('Maximum error for smartphone 1: ',max(err1))
print('Maximum error for smartphone 2: ',max(err2))
print('Minimum error for smartphone 1: ',min(err1))
print('Minimum error for smartphone 2: ',min(err2))
print('Standard deviation for error for smartphone 1: ',np.std(err1))
print('Standard deviation for error for smartphone 2: ',np.std(err2))


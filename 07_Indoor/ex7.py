# -*- coding: utf-8 -*-
"""
@author: Mattia
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import data
user=pd.read_csv('user_db.txt',header=None)
cp=pd.read_csv('control_points_db.txt',header=None)
#cp coordinates
cpx=cp[1]
cpy=cp[2]
#user coordinates
ux=[0]*25
uy=[0]*25
#cp measurements
cp=cp.set_index(0).drop(columns=[13,1,2])
cp=cp.rename(columns={3:1,4:2,5:3,6:4,7:5,8:6,9:7,10:8,11:9,12:10})
#user measurements
user=user.set_index(0).drop(columns=11)
#user coordinates estimation
dist=[0]*36
for i in range (25):
    t=user.iloc[[i]]
    for j in range (36):
        dist[j]=np.dot((user.iloc[[i]]-cp.iloc[j]),np.transpose(user.iloc[[i]]-cp.iloc[j]))
    minpos=dist.index(min(dist))
    ux[i]=cpx[minpos]
    uy[i]=cpy[minpos]
#plot
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(cpx, cpy, '.b', label='Control points')
ax.plot(ux, uy, 'g', label='User trajectory')
plt.xlim([0,20])
plt.ylim([0,10])
plt.title("User trajectory")
leg = ax.legend();
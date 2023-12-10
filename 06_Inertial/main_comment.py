# -*- coding: utf-8 -*-
"""
POSITIONING & LOCATION BASED SERVICES
AA 2021/2022

EX 6: Inertial Navigation
@author: Marianna Alghisi

-------------------------------------------------------------------------------
Guidelines: a vehicle that is moving in the 2D X-Y inertial planar system with 
two accelerometers and a gyroscope mounted on it. According to the provided data
in Body System compute the position of the vehicle in Iertial System

Input data: both .dat files contain:
    - col1: Epoch [s]
    - col2: Acceleration in X [m/s]
    - col3: Acceleration in Y [m/s]
    - col4: Angular velocity [rad/s]
'Inertial_data.dat' --> observations simulated without errors
'Inertial_data_ni.dat' --> observations simulated with errors

Workflow:
    1) Import the data from both datasets (with and without errors)
       HINT: numpy.loadtxt() function can be used to load .dat files as arrays
    2) Use calcTrajectory function to compute the position of the vehicle
       in each epoch
    3) Plot the results 
-------------------------------------------------------------------------------
"""

from re import A
import numpy as np
import matplotlib.pyplot as plt

from function import calcTrajectory

P0 = np.array([[100],[100]])
#1) Import the data
#2) Use the function to compute the trajectory in Inertial RS for both datasets
epoch, a_x, a_y, omegaz= np.loadtxt('Inertial_data.dat',unpack=True)
res_without = calcTrajectory(P0 , a_x, a_y, omegaz, epoch)

epoch, a_x, a_y, omegaz= np.loadtxt('Inertial_data_errors.dat',unpack=True)
res_with = calcTrajectory(P0 , a_x, a_y, omegaz, epoch)

res_without = np.array(res_without)
res_with = np.array(res_with)

#3) Plot the results
fig, ax = plt.subplots (figsize = (10,5))
ax.set(xlabel='x[m]', ylabel='y[m]', title = 'trajectory')
blue, = ax.plot(res_without[:,0], res_without[:,1], '-', color='blue', label = 'No noise')
red, = ax.plot(res_with[:,0], res_with[:,1], '-', color='red', label = 'With noise')


ax.legend(handles=[blue, red])
# ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0)) #scientific notation
plt.savefig('trajectory.png')


# fig, ax = plt.subplots(figsize=(10,6))
# ax.plot(x, y, '-', color='blue')
# plt.savefig('time-height.jpg')
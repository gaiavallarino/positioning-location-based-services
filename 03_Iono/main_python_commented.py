# -*- coding: utf-8 -*-
"""
POSLBS - 3rd Lab: Ionospheric delay computation

@author: Marianna Alghisi
"""

'''
Goals: 
   1) 4 zenithal maps of ionospheric error corrections:
   - Elevation: 90°
   - Latitude: [-80°, 80°, step = 0.5°]
   - Longitude: [-180°, 180°, step = 0.5°]
   - time: [00:00, 06:00, 12:00, 18:00]
   
   2) 2 polar maps of ionospheric error corrections:
    - Observer located in Milan
    - Elevation: [0, 90°, step = 0.5°]
    - Azimuth: [-180°, 180°, step = 0.5°]
    time: [00:00, 12:00]
'''

# Import required libraries
from enum import auto
import math
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from requests import NullHandler
from sqlalchemy import null
from ionoCorrection import calculate_iono_grid,calculate_mi_grid

# Ionospheric correction parameters:
alpha = [7.4506*10**(-9), 1.4901*10**(-8), -5.9605*10**(-8), -1.1921*10**(-7)]
beta = [9.2160*10**(4), 1.3107*10**(5), -6.5536*10**(4), -5.2429*10**(5)]

'''
1) ZENITHAL MAPS
'''
# Initialization of the parameters: define inputs for the Zenithal maps (elevation, azimuth, time, lat, lon)

# Loop on time, latitude and longitude --> compute for each point the Ionospheric delay
# TIP: store latitude, longitude and iono_delay in list objects!

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '06:00', '12:00', '18:00']
GPStimes = [(int(x[0:2])*3600) for x in timePrint] # CHECK LIST COMPREHENSION FOR PYTHON IN GOOGLE IF U DONT UNDERSTAND

lats  = np.arange(-80,80.5,0.5)
longs = np.arange(-180,180.5,0.5)
grid  = np.array(np.meshgrid(lats, longs)).T.reshape(-1,2)

for i in range(4):
    t = timePrint[i]
    IONO = calculate_iono_grid(grid,GPStimes[i],alpha,beta)
    results = pd.DataFrame()
    results['latitude'] = np.transpose(grid)[0]
    results['longitude'] = np.transpose(grid)[1]
    results['iono_delay'] = IONO
    gdf = gpd.GeoDataFrame(results, geometry=gpd.points_from_xy(results.longitude, results.latitude), crs = 3857)
    world = gpd.read_file('world/world.shp')
    fig, ax = plt.subplots (figsize = (15,15))
    world.boundary.plot(ax=ax, color='black')
    ax.set(xlabel='Longitude', ylabel='Latitude', title='Zenithal map of ionospheric delay at '+ str(t))
    gdf.plot(column='iono_delay', ax = ax, marker='o', legend=True)
    fig.savefig(str(t[0:2]) + '.png', dpi=fig.dpi)


'''
2) POLAR MAPS
'''
# Definition of Milano's coordinate
lat_mi = 45 + 28/60 + 38.28/60**2
lon_mi = 9 + 10/60 + 53.4/60**2

# Inizialization of the parameters for the loop: time, elevation, azimuth
elevation  = np.arange(0,90.5,0.5)
azimuth = np.arange(-180,180.5,0.5)
grid_mi  = np.array(np.meshgrid(azimuth, elevation)).T.reshape(-1,2)

timePrint_mi = ['00:00', '12:00']
GPStimes_mi = [(int(x[0:2])*3600) for x in timePrint_mi]

vminimum = math.inf
vmaximum = -math.inf
delays = []
for i in range(0,2):
    delays.append(calculate_mi_grid(lat_mi,lon_mi,GPStimes_mi[i],alpha,beta,grid_mi))
    vminimum = min(vminimum,min(delays[i]))
    vmaximum = max(vmaximum,max(delays[i]))

norm = mpl.colors.Normalize(vmin=vminimum,vmax=vmaximum)

for i in [0,1]:
    t = timePrint_mi[i]
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(projection='polar')
    plt.scatter(np.transpose(grid_mi)[0], np.transpose(grid_mi)[1], c=delays[i], cmap='brg',norm=norm,  alpha=0.75, label=delays)
    ax.set_title('Ionospheric Error Polar Map for Milan Observer time = '+str(t))
    ax.set_rmax(90)
    plt.colorbar(label='Ionospheric Delay')
    fig.savefig(str(t[0:2])+'polar'+'.png', dpi=fig.dpi)

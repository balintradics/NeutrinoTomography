#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 16:52:48 2020

@author: balintradics
"""
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='ortho',lat_0=0,lon_0=0,resolution='l')
#map = Basemap(projection='moll',lon_0=0,resolution='c')
#map = Basemap(projection='mill',lon_0=180)


# draw coastlines, country boundaries, fill continents.
# add an axes with a black background
fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
map.drawcoastlines(linewidth=0.25)
map.shadedrelief(scale=0.5)
#map.drawcountries(linewidth=0.25)
map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
map.fillcontinents(color='coral',lake_color='aqua')

#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
#map.drawmeridians(np.arange(0,360,30))
#map.drawgreatcircle(30, -90, 30, 90)

# make up some data on a regular lat/lon grid.
f = open('/Users/balintradics/Work/Neutrino/EarthViz/flux_map.dat')
tables = np.genfromtxt(f)
f.close()

R_E = 6371. # km
dCell = 300. # km
delta = dCell/R_E 

Lat = [] # Latitude
Lon = [] # Longitude
A = [] # values

Lat = tables[:,0]
Lon = tables[:,1]
A = tables[:,2]
Zmin = min(A)
Zmax = max(A)

# create an array with [x_0, y_0, z_0]
rows = []
for i in range(len(Lat)):
    rows.append([Lon[i], Lat[i], A[i]])
    
mat = np.array(rows)
    
# Create grid on which we interpolate in 2d
lats_new = np.arange(-0.5*np.pi, 0.5*np.pi+delta, delta)
lons_new = np.arange(0, 2*np.pi+delta, delta)
lons, lats = np.meshgrid(lons_new, lats_new)



# Interpolating function feeded from data
#f = interpolate.interp2d(Lon, Lat, A, kind='linear')
# interpolate: 'nearest', 'linear', 'cubic'
Z = interpolate.griddata((mat[:,0], mat[:,1]), mat[:,2], (lons,lats), method='nearest')



# compute native map projection coordinates of lat/lon grid.
x, y = map(lons*180./np.pi, lats*180./np.pi)



p = map.pcolor(x,y,Z, cmap='RdBu_r')
#map.contourf(x,y, Z)

plt.colorbar(p)
plt.clim(0.3, 0.5)
#plt.savefig('line_plot_test.pdf')

plt.show()


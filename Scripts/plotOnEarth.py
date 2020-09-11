from mpl_toolkits.basemap import Basemap
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='ortho',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
# add an axes with a black background
fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
map.drawcoastlines(linewidth=0.25)
map.shadedrelief(scale=0.5)
#map.drawcountries(linewidth=0.25)
#map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))
#map.drawgreatcircle(30, -90, 30, 90)

# make up some data on a regular lat/lon grid.
f = open('/Users/balintradics/Work/Neutrino/EarthViz/flux_map.dat')
tables = np.genfromtxt(f)
f.close()

R_E = 6371. # km
dCell = 500. # km
delta = dCell/R_E 

Lat = [] # Latitude
Lon = [] # Longitude
A = [] # values

Lat = tables[:,0]
Lon = tables[:,1]
A = tables[:,2]/1.0e+16
    
# Create grid on which we interpolate in 2d
lats_new = np.arange(-0.5*np.pi, 0.5*np.pi+delta, delta)
lons_new = np.arange(0, 2*np.pi+delta, delta)
lons, lats = np.meshgrid(lons_new, lats_new)

# Interpolating function feeded from data
f = interpolate.interp2d(Lat, Lon, A, kind='linear')

# Interpolated values
wave = f(lons_new, lats_new)


# compute native map projection coordinates of lat/lon grid.
x, y = map(lons*180./np.pi, lats*180./np.pi)

cs = map.pcolor(x,y,wave, cmap='RdBu_r')

plt.colorbar()
plt.clim(1.5, 4)
#plt.savefig('line_plot_test.pdf')

plt.show()

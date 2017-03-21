# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 18:04:27 2016

@author: imchugh
"""

import os
import numpy as np
import linecache
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
import gdal

def truncate_colormap(cmap, minval = 0.0, maxval = 1.0, n = 100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n = cmap.name, a = minval, b = maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

# Get the data
path = '/media/Data/Dropbox/Data_sites non flux/Site data plots and non-flux/Sites/Whroo/Data/DEM/DEM_1sec'
f = 'whroo_dem.asc'
target = os.path.join(path, f)

# Specify tower location
tower_lat = -36.673215
tower_lon = 145.029247

# Read into array with gdal package and get required metadata
ds = gdal.Open(target)
data = ds.ReadAsArray()
stats = ds.GetGeoTransform()
westernmost_lon = stats[0]
northernmost_lat = stats[3]
cell_size_degrees = stats[1]
cells_lat = data.shape[0]
cells_lon = data.shape[1]
easternmost_lon = westernmost_lon + cell_size_degrees * cells_lon
southernmost_lat = northernmost_lat - cell_size_degrees * cells_lat
centre_lat = (southernmost_lat + northernmost_lat) / 2
centre_lon = (easternmost_lon + westernmost_lon) / 2

# This is same as using gdal simple functions, EXCEPT that gdal uses northernmost
# latitude as base coordinate rather than southernmost
#data = np.loadtxt(target, skiprows=6)
#hdr = [linecache.getline(target, i) for i in range(1,7)]
#values = [float(h.split(" ")[-1].strip()) for h in hdr]
#cols, rows, lx, ly, cell, nd = values

# DEM Raster datum is WGS84; EPSG code 4326
# Is the epsg kwarg in the Basemap class instance specifying a projection?!
map = Basemap(projection = 'merc',
              lat_0 = centre_lat, lon_0 = centre_lon,
              llcrnrlat = southernmost_lat,
              llcrnrlon = westernmost_lon,
              urcrnrlat = northernmost_lat,
              urcrnrlon = easternmost_lon)
#              epsg = 4326)

# Get x and y coordinates of space in map coordinates (m)
x = np.linspace(map.llcrnrx, map.urcrnrx, data.shape[1])
y = np.linspace(map.llcrnry, map.urcrnry, data.shape[0])
if map.llcrnrlat < 0:
    y = y [::-1]

# Create longitude / latitude grid onto which to project the elevation data
xx, yy = np.meshgrid(x, y)

# Do plotting
fig, ax = plt.subplots(1, 1, figsize = (12, 9))
fig.patch.set_facecolor('white')
cmap = plt.get_cmap('cubehelix')
new_cmap = truncate_colormap(cmap, 0.5, 1)
color_min = int(data.min() / 10.0) * 10
color_max = int(data.max() / 10.0 + 1) * 10
colormesh = map.pcolormesh(xx, yy, data, vmin = color_min, 
                           vmax = color_max, cmap=new_cmap)
contour = map.contour(xx, yy, data, colors = 'k')
plt.clabel(contour, inline=True, fmt='%1.0f', fontsize=12, colors='k')
cb = map.colorbar(colormesh, pad = '5%')
cb.set_label('Elevation (m)', labelpad = 10)
clr_ax = cb.ax
text = clr_ax.yaxis.label
font = matplotlib.font_manager.FontProperties(size = 18)
text.set_font_properties(font)

# Mark location of tower
point_x, point_y = map(tower_lon, tower_lat)
map.plot(point_x, point_y, marker = 's', color = 'black', markersize = 8)
plt.text(point_x + 250, point_y - 80, 'Tower', fontsize = 18)

# Create north arrow         
arrow_lon = easternmost_lon - abs(westernmost_lon - easternmost_lon) / 12
arrow_lat = northernmost_lat - abs(northernmost_lat - southernmost_lat) / 15
arrow_x, arrow_y = map(arrow_lon, arrow_lat)
text_x, text_y = (0, -120)
plt.annotate('$N$', xy = (arrow_x, arrow_y), xycoords = 'data',
                xytext = (text_x, text_y), textcoords = 'offset points',
                color='black', fontsize = 24, verticalalignment = 'center',
                horizontalalignment = 'center', 
                arrowprops = dict(shrink = 0.05,
                                  width = 0.5,                                    
                                  color = 'black'))

# Draw scale bar
scale_bar_lon_loc = westernmost_lon + abs(westernmost_lon - easternmost_lon) / 5
scale_bar_lat_loc = southernmost_lat + abs(northernmost_lat - southernmost_lat) / 10
map.drawmapscale(scale_bar_lon_loc, 
                 scale_bar_lat_loc, 
                 tower_lon, tower_lat, 2000, 
                 units = 'm', fontsize = 12, barstyle = 'fancy')

plt.show()

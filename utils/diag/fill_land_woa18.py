import os, sys
import numpy as np
import datetime as dt
import xesmf as xe
from dateutil.relativedelta import relativedelta

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

sys.path.append("../pycommon")
from pps_tools import set_cartopy_gridlines, set_cartopy_colorbar

from regrid_tools import iterative_fill_POP_core

#
# read in the lat & lon of the input WOA18 fields
#
fn_grd_in = 'nudging_utils/woa18_decav_s01_04.nc'
f = Dataset(fn_grd_in)
lat1d_grd_in = f.variables['lat'][:]
lon1d_grd_in = f.variables['lon'][:]
sss_grd_in = np.squeeze(f.variables['s_an'][:])[0,:,:] # use the 1st layer
depth_grd_in = f.variables['depth'][:]
f.close()

lon2d_grd_in, lat2d_grd_in = np.meshgrid(lon1d_grd_in, lat1d_grd_in)


print("lat1d_grd_in: shape, min, max=",lat1d_grd_in.shape, np.min(lat1d_grd_in), np.max(lat1d_grd_in))
print("lon1d_grd_in: shape, min, max=",lon1d_grd_in.shape, np.min(lon1d_grd_in), np.max(lon1d_grd_in))
print("sss_grd_in: shape, min, max=",sss_grd_in.shape, np.min(sss_grd_in), np.max(sss_grd_in))
print("depth:", depth_grd_in)

print("-"*80)
v2d_in = sss_grd_in.copy()
print(type(v2d_in))
print(np.sum(v2d_in.mask), np.sum(v2d_in.mask), np.sum(v2d_in.mask)/(v2d_in.shape[0]*v2d_in.shape[1]))

nlat, nlon = v2d_in.shape

# wk2d is numpy.array, not masked one
wk2d = v2d_in.data  
missing_value = np.float32(1e35) #v2d_in.fill_value  #1e36 #v2d_in.fill_value
print("type(missing_value)=", type(missing_value))
wk2d[v2d_in.mask] = missing_value

print("type(wk2d)=",type(wk2d), np.min(wk2d), np.max(wk2d))

fillmask = v2d_in.mask 
#print("outside: var[48,:]=",v2d_in.data[48,:])
#print("outside: var[48,:]=",v2d_in[48,:])
#print("outside: var[48,:]=",wk2d[48,:])

#
# fill values over land
#
iterative_fill_POP_core(var=wk2d, fillmask=fillmask, missing_value=missing_value, tol=1.e-4, ltripole=True)


#
# final plot check
#


v2d_out = wk2d
# plot the output SSS fields with land masked out
fig = plt.figure()
ax  = fig.add_subplot(111,projection=ccrs.Robinson())
crs = ccrs.PlateCarree()

#surf = ax.contourf(lon2d_grd_in, lat2d_grd_in, v2d_out, levels=np.arange(32,38.5,0.5), extend="both", cmap=plt.cm.jet, transform=crs)
surf = ax.pcolormesh(lon2d_grd_in, lat2d_grd_in, v2d_out, vmin=32, vmax=38.5, edgecolors=None, cmap=plt.cm.jet, transform=crs)
set_cartopy_colorbar(ax,surf,fig,shrink=1)
ax.add_feature(cfeature.COASTLINE)
#ax.add_feature(cfeature.LAND, facecolor=[0.7,0.7,0.7])
set_cartopy_gridlines(ax=ax,crs=crs)
ax.set_title('SSS (unit: psu)',fontsize=14,fontweight='heavy')
plt.show()


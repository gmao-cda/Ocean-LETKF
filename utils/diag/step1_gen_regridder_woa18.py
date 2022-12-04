import numpy as np
import xesmf as xe
from netCDF4 import Dataset

if __name__ == '__main__':
  # read in the output lon, lat, and land/sea mask(sea-1;land-0)
  fn_grd_out = "nudging_utils/flood_salt_restore_PHC2.1440x1080.v20180405.nc"
  f = Dataset(fn_grd_out)
  lat2d_grd_out = f.variables['lat'][:]
  lon2d_grd_out = f.variables['lon'][:]
  f.close()

  #fn_wet = 'nudging_utils/ocean_topog.nc'
  #f = Dataset(fn_wet)
  #wet_grd_out = f.variables['wet'][:]
  #f.close()
  #print("wet: shape, min, max=",wet_grd_out.shape, np.min(wet_grd_out), np.max(wet_grd_out))

  print("lat2d_grd_out: shape, min, max=",lat2d_grd_out.shape, np.min(lat2d_grd_out), np.max(lat2d_grd_out))
  print("lon2d_grd_out: shape, min, max=",lon2d_grd_out.shape, np.min(lon2d_grd_out), np.max(lon2d_grd_out))

  # read in the input lon, lat
  fn_grd_in = 'nudging_utils/woa18_decav_s01_04.nc'
  f = Dataset(fn_grd_in)
  lat1d_grd_in = f.variables['lat'][:]
  lon1d_grd_in = f.variables['lon'][:]
  f.close()

  print("lat1d_grd_in: shape, min, max=",lat1d_grd_in.shape, np.min(lat1d_grd_in), np.max(lat1d_grd_in))
  print("lon1d_grd_in: shape, min, max=",lon1d_grd_in.shape, np.min(lon1d_grd_in), np.max(lon1d_grd_in))

  # generate the regridder 
  grd_in = {"lon": lon1d_grd_in, "lat": lat1d_grd_in}  
  grd_out = {"lon": lon2d_grd_out, "lat": lat2d_grd_out}
  # we are interpolating 0.25 deg SSS to 0.25 tripolar grids
  interp_method = "bilinear"
  periodic = True
  regridder = xe.Regridder(grd_in, grd_out, interp_method, periodic=periodic)
  

  # if regridder is unavailable, take following steps
  writeRegridderIntoFile = True
  fn_regridder_out = "regridInfo_bilinear_720x1440_tripolar_1080x1440_peri.nc"

  if writeRegridderIntoFile:
     print("write regrid info into the file: ", fn_regridder_out)
     regridder.filename = fn_regridder_out
     regridder.to_netcdf()


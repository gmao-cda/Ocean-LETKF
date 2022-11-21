if __name__ == "__main__":
    import os, sys
    import numpy as np
    import h5py
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.util import add_cyclic_point

    from pps_tools import set_cartopy_gridlines, set_cartopy_colorbar


    #
    # load data
    #

    #fndir = '/Users/cda/Documents/work/project/ocean-letkf-develop/MOM6-LETKF-Data/obs/'
    #fname = 'SMAP_L2B_SSS_36950_20220101T005200_R18240_V5.0.h5'
    #fnin= fndir + fname
    fnin = sys.argv[1]
    print(os. path. basename(fnin))
    f = h5py.File(fnin,'r')
    lat2d = f['lat'][:]
    lon2d = f['lon'][:]
    sss = f['smap_sss'][:]
    err_sss = f['smap_sss_uncertainty'][:]
    f.close()

    print("lat2d: shape, min, max=",lat2d.shape, np.min(lat2d), np.max(lat2d))
    print("lon2d: shape, min, max=",lon2d.shape, np.min(lon2d), np.max(lon2d))
    print("sss: shape, min, max=",sss.shape, np.min(sss), np.max(sss))
    print("err_sss: shape, min, max=",err_sss.shape, np.min(err_sss), np.max(err_sss))


    #
    # plot    
    #

    fig  = plt.figure(figsize=(10,8))

    # subplot for SSS 
    ax   = fig.add_subplot(211,projection=ccrs.PlateCarree(central_longitude=180))
    crs = ccrs.PlateCarree()

    surf = ax.contourf(lon2d, lat2d, sss, levels=np.arange(32,38.5,0.5), extend="both", cmap=plt.cm.jet, transform=crs)
    set_cartopy_colorbar(ax,surf,fig,shrink=1)

    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.LAND, facecolor=[0.7,0.7,0.7])
    set_cartopy_gridlines(ax=ax,crs=crs)

    ax.set_title('SSS (unit: psu)',fontsize=14,fontweight='heavy')


    # subplot for uncertainty of SSS
    ax2   = fig.add_subplot(212,projection=ccrs.PlateCarree(central_longitude=180))
    crs = ccrs.PlateCarree()

    #surf = ax2.contourf(lon2d, lat2d, err_sss, levels=np.arange(0,2.0,0.1), extend="max", cmap=plt.cm.jet, transform=crs)
    surf = ax2.scatter(lon2d, lat2d, err_sss, cmap=plt.cm.jet, transform=crs)
    set_cartopy_colorbar(ax2,surf,fig,shrink=1)

    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.LAND, facecolor=[0.7,0.7,0.7])
    set_cartopy_gridlines(ax=ax2,crs=crs)

    ax2.set_title('Error of SSS (unit: psu)',fontsize=14,fontweight='heavy')

    plt.show()
    


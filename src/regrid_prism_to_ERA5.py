"""
regrid_prims_to_ERA5.py
==============================================================
D.T. Milodowski, University of Edinburgh
--------------------------------------------------------------
Code to transform the PRISM precipitation data onto the ERA5
0.125 deg grid.
PRISM data are initially stored as individual daily .bil files
Code structure:
1) For each month...
|   - load ERA5 file as template dataset
|   - subset to ROI
|   - set values to zero
|   2) for each day in month...
|   |   -   regrid to ERA5 lat/lon (WGS84) grid
|   |   -   load regridded file
|   |   -   add daily precipitation value to daily slot in the 
            ERA5 template
|   
|   write dataset to new netcdf file
"""
import numpy as np                  # standard scientific programming library
import xarray as xr                 # geospatial library, very good for netcdfs
import matplotlib.pyplot as plt     # standard plotting library; not necessary for regridding
import os                           # useful package for system commands

# define some paths where data is/will be stored
path2ERA5 = '/exports/csce/datastore/geos/groups/gcel/ECMWF/ERA5/0.125deg_global'
path2PRISM = '/exports/csce/datastore/geos/groups/gcel/jau/PRISM/prismtmp'
path2PRISM_aggregated = '/exports/csce/datastore/geos/groups/gcel/jau/PRISM_0.125deg_CA'

if not os.path.isdir(path2PRISM_aggregated):
    os.mkdir(path2PRISM_aggregated)

# bounding box
N = 42.25
S = 32.25
W = -124.75
E = -108

# target resolution
dx = 0.125
dy = 0.125

# months to do
start = np.datetime64('2014-01')
end = np.datetime64('2023-01')
months = np.arange(start,end)

#---------------------------------------------------------------
# Load in an ERA5 file to get some information about the grid.
# We don't want to load in the whole globe each time, so helps
# to know how to subset the latitude and longitude dimensions.

# Note ERA 5 is "upside down"

# convert date info to text strings
date_string = np.datetime_as_string(start)
mm = date_string[5:]
yy = date_string[:4]

# open netcdf file
era5 = xr.open_dataset('%s/precipitation_daily_mean_%s%s.nc' % (path2ERA5,yy,mm))
lat_mask = np.all((era5.Latitude.values>=S-dy/2.,era5.Latitude.values<=N+dy/2.),axis=0)
lon_mask = np.all((era5.Longitude.values>=W-dx/2.,era5.Longitude.values<=E+dy/2.),axis=0)

lat_idx_0 = era5.lat_dim.values[lat_mask].min()
lat_idx_1 = era5.lat_dim.values[lat_mask].max()
lon_idx_0 = era5.long_dim.values[lon_mask].min()
lon_idx_1 = era5.long_dim.values[lon_mask].max()

# subset to CA
era5 = era5.sel(long_dim=slice(lon_idx_0,lon_idx_1),lat_dim=slice(lat_idx_0,lat_idx_1))

# bounding box for regridding
N_ = era5.Latitude.values.max()+dy/2
S_ = era5.Latitude.values.min()-dy/2
E_ = era5.Longitude.values.max()+dy/2
W_ = era5.Longitude.values.min()-dy/2

#---------------------------------------------------------------
# Loop through the months
for month in months:
    # convert date info to text strings
    date_string = np.datetime_as_string(month)
    mm = date_string[-2:]
    yy = date_string[:4]

    # open up the ERA5 data    
    era5 = xr.open_dataset('%s/precipitation_daily_mean_%s%s.nc' % (path2ERA5,yy,mm)).sel(long_dim=slice(lon_idx_0,lon_idx_1),lat_dim=slice(lat_idx_0,lat_idx_1))
    era5.daily_precip.values*=np.nan # reset to nans so we know if data appears, it is the PRISM data!

    # Loop through the days in the month
    start_of_month = month.astype('datetime64[D]')
    end_of_month = (month+np.timedelta64(1,'M')).astype('datetime64[D]')
    days_in_month = np.arange(start_of_month,end_of_month)
    for ii,day in enumerate(days_in_month):     # enumerate gives the index and the element in the array for each step of the loop
        # get the day as a string
        dd = np.datetime_as_string(day)[-2:]
        # regrid PRISM to 0.125 deg, WGS84, using gdalwarp via command line
        path2bil = '%s/PRISM_ppt_stable_4kmD2_%s%s%s_bil' % (path2PRISM,yy,mm,dd)
        PRISM_file = '%s/PRISM_ppt_stable_4kmD2_%s%s%s_bil.bil' % (path2bil,yy,mm,dd)
        
        # a temp file for the projected raster
        temp_file = '%s/temp_%s%s%s.bil' % (path2PRISM,yy,mm,dd)
        
        # some gdal magic
        os.system('gdalwarp -tr %f %f -te %f %f %f %f -t_srs EPSG:4326 -r average %s %s' % (dx,dy,W_,S_,E_,N_,PRISM_file,temp_file))

        # Open PRISM file for this day
        PRISM = xr.open_rasterio(temp_file).sel(band=1)    # This will throw a deprecation warning, but for now it's ok to use this function
        
        # Set nodata values
        PRISM.values[PRISM.values==-9999] = np.nan # Not really necessary as we will convert back before saving, but useful for visualising if plotting
        
        # fill era5 array
        era5.daily_precip.values[ii] = PRISM.values.copy() # if the arrays have the right dimensions, this should be fine. Safer to use copy() rather than direct assignment

        # remove temp file
        os.system('rm %s' % temp_file)

        # end day loop

    # save monthly output to netcdf
    outfile = '%s/precipitation_daily_mean_%s%s.nc' % (path2PRISM_aggregated,yy,mm)
    era5.to_netcdf(outfile)    
    
    # End month loop
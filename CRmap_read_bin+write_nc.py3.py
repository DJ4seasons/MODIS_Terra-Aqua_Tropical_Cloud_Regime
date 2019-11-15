"""
Read binary 'Cloud Regime RFO Map' file, and write to NetCDF format

Input files:
1. CR_map
data/[aqua or terra]_CRnum_map.MODISc61_b42_DTR_CR10_sub4.from[YYYY]0101_365dx50x360.int16dat
[nt=365 or 366 days, nLat=50 (25S-25N), nLon=360 (-179.5,179.5)], int16 (2 Byte)

-1: Missing
1-9: Regime 1 to 9
101-104: Regime 10; sub-regime 1 to 4

By Daeho Jin
2019.11.14
"""

import sys
import numpy as np
import os.path
from subprocess import call
from datetime import timedelta, date, datetime
from netCDF4 import Dataset, date2num

def bin_file_read2mtx(fname,dtp=np.float32):
    """ Open a binary file, and read data
        fname : file name
        dtp   : data type; np.float32 or np.float64, etc. """

    if not os.path.isfile(fname):
        print("File does not exist:"+fname)
        sys.exit()

    with open(fname,'rb') as fd:
        bin_mat = np.fromfile(file=fd,dtype=dtp)

    return bin_mat

###--- Parameters
lons_new = np.arange(-179.5,180.,1.); nlon=len(lons_new)
lats_new = np.arange(-24.5,25.,1.); nlat=len(lats_new)

dset= 'aqua' #'terra' #'aqua'
###--
indir = "./Data/"
fn_header= '_CRnum_map.MODISc61_b42_DTR_CR10_sub4.'

###--
outdir= indir
for yr in range(2003,2018,1):
    idate=date(yr,1,1)
    edate=date(yr,12,31)
    ndy_yr=(edate-idate).days+1

    fn_crmap = indir+dset+fn_header+'from{}_{}dx{}x{}.int16dat'.format(idate.strftime('%Y%m%d'),ndy_yr,nlat,nlon)
    crmap= bin_file_read2mtx(fn_crmap,dtp=np.int16).reshape(ndy_yr,nlat,nlon)

    print(yr,ndy_yr)

    fn_crmap_out = indir+dset+fn_header+'{}.nc'.format(yr)
    ncfw= Dataset(fn_crmap_out, "w", format="NETCDF4")

    ## Attribution
    ncfw.description = "Regime Number Map of MODIS Deep Tropical Cloud Regimes"
    #ncfw.set_fill_off

    ## Dimensions
    xx=ncfw.createDimension('nlon',nlon)
    yy=ncfw.createDimension('nlat',nlat)
    tt=ncfw.createDimension('time',ndy_yr)

    ## Data
    lonnc = ncfw.createVariable('lon','f4','nlon')
    lonnc[:] = lons_new
    lonnc.units = 'degrees_East'
    lonnc.description = 'Longitudes'

    latnc = ncfw.createVariable('lat','f4','nlat')
    latnc[:] = lats_new
    latnc.units = 'degrees_North'
    latnc.description = 'Latitudes'

    timenc = ncfw.createVariable('time','f8','time')
    timenc.units = 'days since 0001-01-01'
    timenc.calendar = 'gregorian'
    times = [datetime(yr,1,1)+timedelta(days=n) for n in range(ndy_yr)]
    timenc[:] = date2num(times, timenc.units, calendar=timenc.calendar)

    crmapnc = ncfw.createVariable('CRnum','i2',('time','nlat','nlon'))
    crmapnc[:] = crmap
    crmapnc.units = 'n/a'
    crmapnc.description = 'Missing=-1; Regime#1 to #9 = 1-9; Regime#10 = 101-104 for sub-regime #1 to #4'

    ## Close file
    print("{} is written.\n".format(fn_crmap_out))
    ncfw.close()

###--- Check the written file
def open_netcdf(fname):
    if not os.path.isfile(fname):
        print("File does not exist:"+fname)
        sys.exit()

    fid=Dataset(fname,'r')
    print("Open:",fname)
    return fid

nc_f= open_netcdf(fn_crmap_out)
print("\n*** NC Format=",nc_f.data_model)

###--- Attributes
print("\n*** Global Attributes ***")
nc_attrs= nc_f.ncattrs()
for nc_attr in nc_attrs:
    print('   {}: {}'.format(nc_attr,nc_f.getncattr(nc_attr)))

print("\n*** Dimensions ***")
for nc_dim in nc_f.dimensions:
#    print('   Name: {}'.format(nc_dim))
    print('   {}'.format(str(nc_f.dimensions[nc_dim]).split(':')[1]))

print("\n*** Variables ***")
for var in nc_f.variables:
    print(nc_f.variables[var])

nc_f.close()

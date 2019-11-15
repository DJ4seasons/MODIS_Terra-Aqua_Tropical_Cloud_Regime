"""
Read binary centroid file, and write to NetCDF format

Input files:
1. Centroids
data/Centroid.MODIS_T+A_b42_DTR_CR10.dpdat
[k=10, nCTP=7, nTAU=6], float64 (8 Byte)

2. Sub-centroids for Regime#10 (last one)
data/Centroid.MODIS_T+A_b42_DTR_CR10.subCR4.dpdat
[k=4, nCTP=7, nTAU=6], float64 (8 Byte)

By Daeho Jin
2019.11.14
"""

import sys
import numpy as np
import os.path
from subprocess import call
#from datetime import timedelta, date, datetime
from netCDF4 import Dataset #, date2num

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
km=10
sub_km=4
ntau=6; nctp=7

tau_boundaries=[0,1.3,3.6,9.4,23,60,150]
ctp_boundaries=[0,180,310,440,560,680,800,1100]

###--
indir = "./Data/"
fn_ctd = 'Centroid.MODIS_T+A_b42_DTR_CR{:d}.dpdat'.format(km)
fn_ctdsub = 'Centroid.MODIS_T+A_b42_DTR_CR{}.subCR{}.dpdat'.format(km,sub_km)

ctd= bin_file_read2mtx(indir+fn_ctd,dtp=np.float64).reshape(km,nctp,ntau)
subctd= bin_file_read2mtx(indir+fn_ctdsub,dtp=np.float64).reshape(sub_km,nctp,ntau)

#print(ctd.sum(axis=1),subctd.sum(axis=1))
print(ctd[0,0,:])
print(ctd[0,-1,:])

outdir = indir
### Ready to Create a NC file
outfn=outdir+'Centroid.MODIS_T+A_b42_DTR_CR{}+subCR{}.nc'.format(km,sub_km)

ncfw= Dataset(outfn, "w", format="NETCDF4")

## Attribution
ncfw.description = "Centroid of MODIS Deep Tropical Cloud Regimes (Terra & Aqua)"

## Dimensions
xx=ncfw.createDimension('nTAU',ntau)
xxb=ncfw.createDimension('n_TAU_Boundary',ntau+1)
yy=ncfw.createDimension('nCTP',nctp)
yyb=ncfw.createDimension('n_CTP_Boundary',nctp+1)
kk1=ncfw.createDimension('nk',km)
kk2=ncfw.createDimension('nk_sub',sub_km)

## Data
tau_b = ncfw.createVariable('TAU_Bin_Boundary','f4','n_TAU_Boundary')
tau_b[:]=tau_boundaries
tau_b.units = 'n/a'
tau_b.description = 'Cloud Optical Depth Boundaries'

ctp_b = ncfw.createVariable('CTP_Bin_Boundary','f4','n_CTP_Boundary')
ctp_b[:]=ctp_boundaries
ctp_b.units = 'hPa'
ctp_b.description = 'Cloud Top Pressure Boundaries'

ctdnc = ncfw.createVariable('CTD','f8',('nk','nCTP','nTAU'))
ctdnc[:] = ctd
ctdnc.units = 'Cloud Fraction'
ctdnc.description = 'Centroids'

subctdnc = ncfw.createVariable('CTD10_sub','f8',('nk_sub','nCTP','nTAU'))
subctdnc[:] = subctd
subctdnc.units = 'Cloud Fraction'
subctdnc.description = 'Sub Centroids of Cloud Regime#10'

## Close file
print("{} is written.\n".format(outfn))
ncfw.close()

###--- Check the written file
def open_netcdf(fname):
    if not os.path.isfile(fname):
        print("File does not exist:"+fname)
        sys.exit()

    fid=Dataset(fname,'r')
    print("Open:",fname)
    return fid

nc_f= open_netcdf(outfn)
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

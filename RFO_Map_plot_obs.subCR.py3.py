"""
#
# Read the CRnum file, and display RFO map
# for 42 histogram bins clustering results
#
# For Tropical Regimes, Sub-CRs
# By Daeho Jin
# 2019.11.14
"""

import numpy as np
import sys
import os.path
#from subprocess import call
from datetime import timedelta, date

#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.colors as cls
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, MultipleLocator
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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

def map_common(ax,label_idx=[True,True,False,True]):
    ax.set_extent([0,359.9,-15.1,15.1],ccrs.PlateCarree())

    ax.coastlines(color='silver',linewidth=1.)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.6, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = label_idx[2]
    gl.ylabels_left = label_idx[0]
    gl.ylabels_right = label_idx[1]
    if ii<nk-1:
        gl.xlabels_bottom = label_idx[3]

    gl.xlocator = FixedLocator(range(-120,361,60)) #[0,60,180,240,360]) #np.arange(-180,181,60))
    gl.ylocator = MultipleLocator(15)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 11, 'color': 'k'}
    gl.ylabel_style = {'size': 11, 'color': 'k'}

    ax.set_aspect('auto')

def add_colorbar(ax):
    ### Add colorbar
    ###- Get position from previous subplot
    pos1=ax.get_position().bounds  ##<= (left,bottom,width,height)
    cb_ax = fig.add_axes([0.1,pos1[1]-pos1[3]/1.8,0.8,pos1[3]/5.])

    tt=np.arange(0,76,15)
    tt2=[str(x)+'%' for x in tt]
    cb = fig.colorbar(cs,cax=cb_ax,orientation='horizontal',ticks=tt,extend='max')
    cb.ax.set_xticklabels(tt2,size=12)

###-- Parameters and defalut values
###-------------------------------------

nk=10; sub_nk=4; nelem=42
date_range=[date(2003,1,1),date(2017,12,31)]
nlat0,nlon0 = 50,360

lons_new = np.arange(-179.5,180.,1.); nlon=len(lons_new)
lats_new = np.arange(-14.5,15.,1.); nlat=len(lats_new)
lons2d, lats2d = np.meshgrid(lons_new,lats_new)
latidx=[10,40]

dir1 = './Data/'

###
###---- Read Centroid
ctdfnm = 'Centroid.MODIS_T+A_b42_DTR_CR{}.subCR{}'.format(nk,sub_nk)
fnm=dir1+ctdfnm+'.f64dat'
obsctd = bin_file_read2mtx(fnm,dtp=np.float64).reshape([sub_nk,nelem])
obscf = np.sum(obsctd,axis=1)*100.

np.set_printoptions(precision=3,suppress=True)
print(obscf)

###
###---- Read CRnum file
fnm1='terra_CRnum_map.MODISc61_b42_DTR_CR{}_sub{}.'.format(nk,sub_nk) #'.5844dx180x360.int16dat'
fnm2='aqua_CRnum_map.MODISc61_b42_DTR_CR{}_sub{}.'.format(nk,sub_nk) #'.5844dx180x360.int16dat'

map0t=np.zeros([sub_nk,nlat,nlon],dtype=float)
map0a=np.zeros([sub_nk,nlat,nlon],dtype=float)

iyr,eyr=date_range[0].year,date_range[1].year
totdays=0
for yy in range(iyr,eyr+1,1):
    idate = date(yy,1,1)
    ndays = (date(yy,12,31)-idate).days+1
    fnm_tail = "from{}_{}dx{}x{}.int16dat".format(idate.strftime('%Y%m%d'),ndays,nlat0,nlon0)

    crnum_t = bin_file_read2mtx(dir1+fnm1+fnm_tail,dtp=np.int16).reshape([ndays,nlat0,nlon0])[:,latidx[0]:latidx[1],:]
    crnum_a = bin_file_read2mtx(dir1+fnm2+fnm_tail,dtp=np.int16).reshape([ndays,nlat0,nlon0])[:,latidx[0]:latidx[1],:]
    for i in range(1,sub_nk+1,1):
        idx1 = crnum_t==i+100
        idx2 = crnum_a==i+100

        map0t[i-1,:,:]+= idx1.sum(axis=0)
        map0a[i-1,:,:]+= idx2.sum(axis=0)

    totdays+=ndays

map0t = map0t/totdays*100.
map0a = map0a/totdays*100.

grfoobst=map0t.reshape([-1,nlat*nlon]).mean(axis=1)
grfoobsa=map0a.reshape([-1,nlat*nlon]).mean(axis=1)
print(grfoobst)
print(grfoobsa)

map0= (map0t+map0a)/2.0
grfoobsall= map0.reshape([-1,nlat*nlon]).mean(axis=1)

###-------------------------------------

###-- Plotting basics
#fig, axs = plt.subplots(3,3)  ## (ny,nx)
fig = plt.figure()
fig.set_size_inches(7,13)    ## (lx,ly)
plt.suptitle("MODIS C6.1 T+A TCR10 Sub-CR RFO Map",fontsize=18,y=1.0)
lf=0.05;rf=0.95
bf=0.03;tf=0.95
fig.subplots_adjust(hspace=0.45,wspace=0.04,left=lf,right=rf,top=tf,bottom=bf)

cm = plt.cm.get_cmap('CMRmap_r',100) #'CMRmap_r' 'YlOrBr' 'Accent' 'afmhot_r'
cmnew = cm(np.arange(100))
cmnew = cmnew[1:96,:] #print cmnew[0,:],cmnew[-1,:]
newcm = cls.LinearSegmentedColormap.from_list("newCMR",cmnew)
newcm.set_under("white")
props = dict(cmap=newcm,alpha=0.9,transform=ccrs.PlateCarree()) #,projection=ccrs.PlateCarree())
abc='abcd'
for ii in range(sub_nk):
    ax=fig.add_subplot(nk,1,ii+1,projection=ccrs.PlateCarree(central_longitude=180))
    # add a title.
    subtit='CR10{} [CF={:.1f}%]: RFO={:.1f}% [T{:.1f}%, A{:.1f}%]'.format(abc[ii],obscf[ii],grfoobsall[ii],grfoobst[ii],grfoobsa[ii])
    print(subtit)
    ax.set_title(subtit,x=0.0,ha='left',fontsize=13,stretch='semi-condensed')

    cs=ax.contourf(lons_new,lats_new,map0[ii,:,:],np.linspace(0,75,101),**props)

    if ii<nk-1:
        label_idx=[True,True,False,False]
    else:
        label_idx=[True,True,False,True]
    map_common(ax,label_idx)

add_colorbar(ax)

###-----------------------------------
### Show or Save
outdir = "./Pics/"
fnout = "RFO."+ctdfnm+".png"

#plt.show()
plt.savefig(outdir+fnout,bbox_inches='tight',dpi=175)
#plt.savefig(outdir+fnout,dpi=160)


#if os.path.isfile(outdir+fnout) and not os.path.isfile('./Pics/'+fnout):
#    call(["ln","-s",outdir+fnout,"./Pics/"])

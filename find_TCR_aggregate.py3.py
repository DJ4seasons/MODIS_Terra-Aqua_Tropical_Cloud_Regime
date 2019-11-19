"""
#
# Read the CRnum file, and identify TCR1 or TCR123 aggregates
# * 42-bin clustering results *
# * MODIS Observation TR10+sub4*
#
# Output: Text File
# [Date, Center of TCR1 (lat [deg N], lon [deg E]), Count of TCR1 (, TCR2, TCR3 if applicapable)]
#
# Daeho Jin, 2019.11.14
#
"""

import numpy as np
from subprocess import call
import os.path
import sys
from datetime import timedelta, date

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

def find_neighbor(ixy,clst,amap,amap2):
    """
    Find neighbor for 4-direction; Recursive
    ixy: tuple, (iy,ix)
    amap: an index map of a given day, T/F
    amap2: an index map of next day, T/F
    clist: xys for one group
    """

    global nlon
    nxys=neighbor_info(ixy)
    for y1,x1 in nxys:
        if x1<0:
            if amap2[y1,x1+nlon]:
                clst.append((y1,x1))
                amap2[y1,x1]=False
                #amap,amap2=
                find_neighbor((y1,x1),clst,amap,amap2)
        else:
            if amap[y1,x1]:
                clst.append((y1,x1))
                amap[y1,x1]=False
                #amap,amap2=
                find_neighbor((y1,x1),clst,amap,amap2)

    return #amap,amap2

def neighbor_info(lxy):
    """
    Input: coordinate of 1 point, (ly,lx)
    Output: a list of coordinate which is neighboring to input point.
          : Only check 4-directions
          : negative lx is allowed (will be used for next day map), but NOT for eastern side, lx can not >=nlon
    """

    global nlon, nlat
    ly,lx=lxy
    nlist=[]
    if ly>0:
        nlist.append((ly-1,lx))
    if ly<nlat-1:
        nlist.append((ly+1,lx))

    nlist.append((ly,lx-1))
    if lx<nlon-1:
        nlist.append((ly,lx+1))

    return nlist

def find_lon_center(array1d):
    global nlon
    xlon=array1d.mean()
    if xlon>nlon:
        xlon-=nlon

    return xlon

def get_cr1_center(cr1s):
    cr1s=np.asarray(cr1s)
    cr1center=[get_latdeg(cr1s[:,0].mean()),
               get_londeg(find_lon_center(cr1s[:,1]))]
    return cr1center

###-------------------------------------

###!!! Parameters and Common parts
options= [('123',5),('1',0)] ## Target TCR(s) and minimum size
opt_idx= 0

tgt_cr = options[opt_idx][0]
sz_crt = options[opt_idx][1]

ncl = 10; subncl=4

###!! Target Dates
tgt_dates= (date(2003,1,1),date(2017,12,31))
tot_days=(tgt_dates[1]-tgt_dates[0]).days+1

lons_new = np.arange(-179.5,180.,1.); nlon=len(lons_new)
lats_new = np.arange(-24.5,25.,1.); nlat=len(lats_new)
#lons2d, lats2d = np.meshgrid(lons_new,lats_new)
get_londeg = lambda x: -179.5+x if x>=180 else 180.5+x
get_latdeg = lambda y: -24.5+y

###!! Get Obs RFO map
dir_obs = './Data/'
fnm_a_head=dir_obs+'aqua_CRnum_map.MODISc61_b42_DTR_CR{}_sub{}.'.format(ncl, subncl) #from{}_{}dx{}x{}.int16dat'.format(t1date.strftime('%Y%m%d'),ndy,nlat,nlon)
fnm_t_head=dir_obs+'terra_CRnum_map.MODISc61_b42_DTR_CR{}_sub{}.'.format(ncl, subncl)

iyr,eyr=tgt_dates[0].year,tgt_dates[1].year
totdays=0
for yy in range(iyr,eyr+1,1):
    idate = date(yy,1,1)
    ndays = (date(yy,12,31)-idate).days+1
    fnm_tail = "from{}_{}dx{}x{}.int16dat".format(idate.strftime('%Y%m%d'),ndays,nlat,nlon)

    crnum_a = bin_file_read2mtx(fnm_a_head+fnm_tail,dtp=np.int16).reshape([ndays,nlat,nlon])
    crnum_t = bin_file_read2mtx(fnm_t_head+fnm_tail,dtp=np.int16).reshape([ndays,nlat,nlon])

    if yy==iyr and idate!=tgt_dates[0]:
        idy= (tgt_dates[0]-idate).days
        crnum_a = crnum_a[idy:,:,:]
        crnum_t = crnum_t[idy:,:,:]

    if yy==eyr and date(yy,12,31)!=tgt_dates[1]:
        edy= (date(yy,12,31)-tgt_dates[1]).days
        crnum_a = crnum_a[:edy,:,:]
        crnum_t = crnum_t[:edy,:,:]

    idx_ms= crnum_a==-1  ## Where missings are in Aqua
    crnum_a[idx_ms] = crnum_t[idx_ms]  ## Filled by Terra

    if yy==iyr:
        crmap = np.copy(crnum_a)
    else:
        crmap = np.concatenate((crmap,crnum_a),axis=0)

print(crmap.shape)
## Add one more time step with full-missings
crmap=np.concatenate((crmap,np.full([1,nlat,nlon],-1)),axis=0)
print(crmap.shape)

###!! Loop for each day
outdir='./Data/'
outnm=outdir+'horizontal_cluster_size+center.MODIS_Aqua+T_TR10.b42_TCR{}.gt{}.{}-{}_25S-25N.txt'.format(tgt_cr,sz_crt,tgt_dates[0].strftime('%Y%m%d'),tgt_dates[1].strftime('%Y%m%d'))
f=open(outnm,'w')
header="date, TCR1 Center(lat_deg,lon_deg), Size of TCR1 (,2 ,3)"
f.write(header+"\n")

for k in range(tot_days):
    tdate=tgt_dates[0]+timedelta(days=k)  ## Target date

    ## This day's data
    if k==0:
        crmap1 = crmap[k,:,:]

        cr1_amap = crmap1==1
        if tgt_cr=='1':
            amap = cr1_amap
        elif tgt_cr=='123':
            amap = np.logical_or.reduce((cr1_amap,crmap1==2,crmap1==3))
        else:
            sys.exit('Other tgt_cr is not working.')
    else:
        crmap1 = crmap2[:]
        cr1_amap = cr1_amap2[:]
        amap = amap2[:]

    ## Next day's data
    crmap2 = crmap[k+1,:,:]
    cr1_amap2 = crmap2==1
    if tgt_cr=='1':
        amap2 = cr1_amap2
    elif tgt_cr=='123':
        amap2 = np.logical_or.reduce((cr1_amap2,crmap2==2,crmap2==3))

    ## Loop untill all aggregates are identified
    while cr1_amap.sum()>0:
        aggr=[]
        where_true= np.where(cr1_amap==True)
        init_xy = [where_true[0][0],where_true[1][0]]
        aggr.append(init_xy)

        cr1_amap[init_xy[0],init_xy[1]]=False

        amap[init_xy[0],init_xy[1]]=False
        find_neighbor(init_xy,aggr,amap,amap2)

        idx_false = amap==False; cr1_amap[idx_false]=False
        idx_false = amap2==False; cr1_amap2[idx_false]=False

        nn=len(aggr)
        if nn>sz_crt:
            ncount=np.zeros([3,],dtype=int)
            cr1s=[]
            for gcell in aggr:
                iy,ix=gcell
                if ix>=0:
                    cridx = crmap1[iy,ix]
                else:
                    cridx = crmap2[iy,ix]
                ncount[cridx-1]+=1
                if cridx==1:
                    cr1s.append(gcell)
            cr1center= get_cr1_center(cr1s)
            if nn!=ncount.sum():
                print('Unexpected Size',tdate,nn,ncount)
                sys.exit()

            outtxt0="{},".format(tdate.strftime('%Y-%m-%d'))
            outtxt1="{:.2f},{:.2f},".format(*cr1center)
            outtxt2="{:d},{:d},{:d}".format(*ncount)
            f.write(outtxt0+outtxt1+outtxt2+"\n")

    if (k+1)%365==0:
        print(k+1,cr1center,ncount) #,nn,max(tmp)

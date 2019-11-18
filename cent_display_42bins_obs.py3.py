"""
#
# Read the centroid file, and display centroids
# for 42 histogram bins clustering results
#
# By Daeho Jin
# 2019.11.14
"""

import numpy as np
import sys
import os.path
#from subprocess import call

#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.colors as cls
import matplotlib.pyplot as plt

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

def cent_show(ax1,ctd,props):
    global ntau,nctp

    pic1=ax1.imshow(ctd,interpolation='nearest',**props)

    ### Axis Control
    ax1.set_xlim(-0.5,ntau-0.5)
    ax1.set_ylim(-0.5,nctp-0.5)
    ax1.set_xticks(np.arange(ntau+1)-0.5)
    ax1.set_yticks(np.arange(nctp+1)-0.5)

    for j in range(7):
        for i in range(6):
            if abs(ctd[j,i])>4.5:
                #ax1.annotate(str(ctd[j,i]),xy=(ix,iy))
                ax1.annotate("{:.0f}".format(ctd[j,i]),xy=(i,j),ha='center',va='center',stretch='semi-condensed',fontsize=11)
    return pic1

def cent_show_common(ax1,subtit,xlabs,ylabs,ylpos='left'):

    print(subtit)
    ax1.set_title(subtit,x=0.0,ha='left',fontsize=14,stretch='semi-condensed')

    ax1.set_xticklabels(xlabs)
    if ylpos=='none':
        ax1.set_yticklabels([])
    else:
        ax1.set_yticklabels(ylabs)

    if ylpos=='right':
        ax1.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

    ax1.yaxis.set_ticks_position('both')

    ### Draw Guide Line
    ax1.axvline(x=1.5,linewidth=0.7,color='k',linestyle=':')
    ax1.axvline(x=3.5,linewidth=0.7,color='k',linestyle=':')
    ax1.axhline(y=1.5,linewidth=0.7,color='k',linestyle=':')
    ax1.axhline(y=3.5,linewidth=0.7,color='k',linestyle=':')

    ### Ticks
    ax1.tick_params(axis='both',which='major',labelsize=11)
    ax1.set_aspect('auto')
    return

def add_colorbar(ax1,pic1):
    ### Add colorbar
    tt=[0.1,0.3,1,3,10,30]
    tt2=[str(x)+'%' for x in tt]

    ###- Get position from previous subplot
    pos1=ax1.get_position().bounds  ##<= (left,bottom,width,height)
    new_left=pos1[0]+pos1[2]*1.15
    cb_ax = fig.add_axes([new_left,pos1[1]-pos1[3]*0.05,pos1[2]*0.1,pos1[3]*1.1])

    cb = fig.colorbar(pic1,cax=cb_ax,orientation='vertical',ticks=tt,extend='both')
    cb.ax.set_yticklabels(tt2,size=12,stretch='semi-condensed')
    cb.ax.minorticks_off()


###-- Parameters and defalut values
nk=10
ntau,nctp=6,7
nelem=ntau*nctp

tau_boundaries=[0,1.3,3.6,9.4,23,60,150]
ctp_boundaries=[1100,800,680,560,440,310,180,0]

###-------------------------------------
dir1 = './Data/'
ctdfnm = 'Centroid.MODIS_T+A_b42_DTR_CR{}'.format(nk)

fnm=dir1+ctdfnm+'.f64dat'
#subfnm=dir1+ctdfnm+'_subCR2.f64dat'

###
###---- Read Centroid
obsctd = bin_file_read2mtx(fnm,dtp=np.float64).reshape([nk,nelem])
obscf = np.sum(obsctd,axis=1)*100.
np.set_printoptions(precision=3,suppress=True)
print(obscf)

###-------------------------------------

###-- Plotting basics
#fig, axs = plt.subplots(4,3)  ## (ny,nx)
fig = plt.figure()
fig.set_size_inches(7.5,11)    ## (lx,ly)
plt.suptitle("MODIS C6 T+A Tropical Cloud Regime",fontsize=20,y=0.98)

## Panel Size
iix=0.05; gapx=0.045; npnx=3
lx=(0.95-iix-gapx*(npnx-1))/float(npnx)
iiy=0.92; gapy=0.055; npny=4
ly=(iiy-0.05-gapy*(npny-1))/float(npny)

ix=iix; iy=iiy

## Color Control
cm = plt.cm.get_cmap('jet',512)
cmnew = cm(np.arange(512))
cmnew = cmnew[72:,:]    #print cmnew[0,:],cmnew[-1,:]
newcm = cls.LinearSegmentedColormap.from_list("newJET",cmnew)
newcm.set_under('white')

props = dict(norm=cls.LogNorm(vmin=0.1,vmax=30),cmap=newcm,alpha=0.9)

for ii in range(nk):
    vv=np.copy(obsctd[ii,:].reshape([nctp,ntau]))*100.
    vv=vv[::-1,:]  # Low CTP going to Top
    ax=fig.add_axes([ix,iy-ly,lx,ly])
    pic1=cent_show(ax,vv,props)

    subtit= "CR{}. CF={:.1f}%".format(ii+1,obscf[ii])
    if ii%3==0:
        ylpos='left'
    elif ii%3==2:
        ylpos='right'
    else:
        ylpos='none'
    cent_show_common(ax,subtit,tau_boundaries,ctp_boundaries,ylpos)

    if ii==int((nk-1)/3)*3:
        ax.set_xlabel('Optical Thickness',fontsize=13)
        ax.set_ylabel('Pressure (hPa)',fontsize=13,labelpad=0)
    if ii==nk-1:
        add_colorbar(ax,pic1)

    ix=ix+lx+gapx
    if ix+lx>1.:
        ix=iix; iy=iy-ly-gapy

###-----------------------------------
### Show or Save
outdir = './Pics/'
fnout = ctdfnm+".png"

#plt.show()
plt.savefig(outdir+fnout,bbox_inches='tight',dpi=175)
#plt.savefig(outdir+fnout,dpi=160)


#if os.path.isfile(outdir+fnout) and not os.path.isfile('./Pics/'+fnout):
#    call(["ln","-s",outdir+fnout,"./Pics/"])

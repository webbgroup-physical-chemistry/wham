#! /anaconda/bin/python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../include')
from mcGenerate_Trajectory import make_traj
#from mcTest_hdf5 import read_hdf5, plot_2d
from read_hdf5 import Read_H5, plot_2d


kb = .0083144621 # kJ/(mol*K)
T = 300 # K
beta = 1/(kb*T)
fc = 50
whamstep = 5
dphi = 0
window_size = 30 # degrees, 12 windows
multiplier = 1 
frames = multiplier*10000
potential = "rb"
whamscript = "../../../bin/wham"
ws = 0
we = 360+ws
datdir = os.path.abspath(os.curdir)
nconv = multiplier*10*100 # check convergence every [nconv] frames
traj = make_traj(kb,T,fc,whamstep,dphi,window_size,frames,potential,ws,we,datdir,stype='mc')

if not os.path.isfile(whamscript) :
    print "ERROR! Cannot find wham binary."
    print "Usage: ./configure && make && make intall"
    sys.exit()

"""
ana1D[0] = angle
ana1D[1] = probability
ana1D[2] = pmf

ana2D[0] = x angle
ana2D[1] = y angle
ana2D[2] = probability
ana2d[3] = pmf
"""
ana1D, ana2D = traj.unbiased_p_distr()
#fl1D, fl2D = traj.get_samples()
fl1D,fl2D = "1D_list.file","2D_list.file"
# Do 1D analysis
outname = "1d.h5"
cmd = "rm %s ; time %s -f %s -o %s -w %s -d 1 -c %i"%(outname, whamscript, fl1D, outname, whamstep, nconv)
os.system(cmd)
g = Read_H5(outname,1)
f,ax = plt.subplots(2)
g.conv_plot(ax[0],ana1D[1])
dat = g.get_dataset('/Ensemble/Probability')
ax[1].plot(dat[0],dat[1],label='WHAM')
ax[1].plot(ana1D[0],ana1D[1],label='analytic')
ax[1].legend(loc='best',ncol=6)
ax[1].set_xlim(-180,175)
ax[1].set_xlabel('Angle (degrees)')
ax[1].set_ylabel('Probability')
ax[1].set_xticks(range(-120,121,60))
f.subplots_adjust(hspace=0.4)
plt.savefig('%s/new_1D.pdf'%datdir,format='pdf')

# Do 2D analysis
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
outname = "2d.h5"
cmd = "rm %s ; time %s -f %s -o %s -w %s -d 2 -c %i "%(outname, whamscript, fl2D, outname, whamstep, nconv)
os.system(cmd)
h = Read_H5(outname,2)
f = plt.figure()
convax = plt.subplot(211)
h.conv_plot(convax,ana2D[2])
whamax = plt.subplot(234)
x,y,z = h.get_dataset('/Ensemble/Probability')
plot_2d(x,y,z,whamax,logscale=True,cb=True,zmin_avg=False)
whamax.set_title('WHAM\nProbability')
whamax.set_ylabel('DoF 2')
anaax = plt.subplot(235)
plot_2d(ana2D[0],ana2D[1],ana2D[2],anaax,logscale=True,cb=True,zmin_avg=False)
anaax.set_title('Analytic\nProbability')
anaax.set_xlabel('DoF 1')
f.subplots_adjust(hspace=0.4,wspace=0.3)
deltaax = plt.subplot(236)
plot_2d(x,y,abs(np.array(ana2D[2]-z)),deltaax,logscale=True,cb=True,zmin_avg=False)
deltaax.set_title('|Difference|')
plt.savefig('%s/new_2D.pdf'%datdir,format='pdf')

"""
cmd = "rm -rvf dat_* 1D_list.file 2D_list.file"
os.system(cmd)
"""

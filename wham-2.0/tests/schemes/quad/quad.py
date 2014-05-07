#! /anaconda/bin/python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../include')
from mcGenerate_Trajectory import make_traj
from mcTest_hdf5 import read_hdf5, plot_2d

kb = .0083144621 # kJ/(mol*K)
T = 300 # K
beta = 1/(kb*T)
fc = 50
whamstep = 5
dphi = 0
window_size = 30 # degrees, 12 windows
frames = 10000
potential = "rb"
whamscript = "../../../bin/wham"
ws = 0
we = 360+ws
datdir = os.path.abspath(os.curdir)
nconv = 100 # check convergence every [nconv] frames
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
fl1D, fl2D = traj.get_samples()
#fl1D,fl2D = "1D_list.file","2D_list.file"
# Do 1D analysis
outname = "1d.h5"
cmd = "rm %s ; %s -f %s -o %s -w %s -d 1 -c %i"%(outname, whamscript, fl1D, outname, whamstep, nconv)
os.system(cmd)
g = read_hdf5(outname,nwindows=360/window_size,ndof=1)
f,ax = plt.subplots(2)
g.converge(ana1D[1],ax[0])
g.plot_prob_1d(ax[1])
ax[1].plot(ana1D[0],ana1D[1],label='analytic')
ax[1].legend(loc='best',ncol=6)
f.subplots_adjust(hspace=0.4)
plt.savefig('%s/1D.pdf'%datdir,format='pdf')

# Do 2D analysis
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)
outname = "2d.h5"
cmd = "rm %s ; %s -f %s -o %s -w %s -d 2 -c %i "%(outname, whamscript, fl2D, outname, whamstep, nconv)
os.system(cmd)
h = read_hdf5(outname,nwindows=360/window_size,ndof=2)#,doconv = False)
f = plt.figure()
convax = plt.subplot(211)
h.converge(ana2D[2],convax)
whamax = plt.subplot(234)
h.plot_prob_2d(whamax,max(ana2D[2]))
whamax.set_title('WHAM\nProbability')
whamax.set_ylabel('DoF 2')
anaax = plt.subplot(235)
plot_2d(ana2D[0],ana2D[1],ana2D[2],anaax)
anaax.set_title('Analytic\nProbability')
anaax.set_xlabel('DoF 1')
f.subplots_adjust(hspace=0.4,wspace=0.3)
deltaax = plt.subplot(236)
h.delta_prob_2d(deltaax,ana2D[2])
deltaax.set_title('|Difference|')
plt.savefig('%s/2D.pdf'%datdir,format='pdf')

cmd = "rm -rvf dat_* 1D_list.file 2D_list.file"
os.system(cmd)

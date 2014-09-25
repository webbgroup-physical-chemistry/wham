#! /usr/bin/env python

"""
    Python equivalent to bw_wham code.  Requires h5py
"""

import sys
import h5py as h5
import numpy as np

try :
    probfile = sys.argv[1]
    filelist = sys.argv[2]
except :
    print "USAGE: %s <h5file> <file list>"%sys.argv[0]
    sys.exit()
    
try :
    f = h5.File(probfile,'r')
except :
    print "Cannot open %s, exiting" %probfile
    sys.exit()
try :
    file = open(filelist,'r')
    filelines = file.readlines()
    file.close()
except :
    print "Cannot open %s, exiting" %filelist
    sys.exit()

# Get the probability of being in each bin
probdat = f['/Ensemble/Probability']
x1,x2,prob = probdat[:,0],probdat[:,1],probdat[:,2]
ntraj = len(f['/Trajectories'].values())

# Read the file list
n = 0
trajectories = {}
xvgs = {}
for line in filelines :
    if not (line.startswith("#") or line.startswith(";") or line.startswith("@")) and line.split() :
        l = line.split()
        trajectories[n] = "/Trajectories/Traj-%i/Bin"%int(l[0])
        xvgs[n] = l[1]
        n += 1

# Make sure the number of xvg files matches the number of trajectories
if n != ntraj :
    print "Error!  The number of files in %s does not match the number of trajectories in %s... Exiting" %(filelist,probfile)
    sys.exit()

# Get the bin assignments for each frame
bins = {}
for i in range(ntraj) :
    bins[i] = f[trajectories[i]][:]

# Prepare to count bins!
nbins = len(prob)
counts = np.zeros(shape=nbins)
# Read each xvg in the file list
ndat = 0
res = {}
var = {}
for i in range(ntraj) :
    print "Reading %s" %xvgs[i]
    try :
        xvg = open(xvgs[i],'r')
        xvglines = xvg.readlines()
        xvg.close()
    except :
        print "Error opening %s (listed in %s), exiting..."%(xvgs[i],filelist)
        sys.exit()
    n = 0
    for line in xvglines :
        if not (line.startswith("#") or line.startswith(";") or line.startswith("@")) and line.split() :
            l = line.split()
            ncol = len(l) - 1
            if ncol > ndat :
                ndat = ncol
            elif ncol < ndat :
                print "Error!  Row %i of %s has %i values, when the previous lines indicated it should have %i values!"%(n,xvgs[i],ncol,ndat)
                sys.exit()
            frame = int(l[0])
            #            print i, frame, bins[i][frame], prob[bins[i][frame]]
            counts[bins[i][frame]] += 1.0
            for j in range(1,1+ncol) :
                try :
                    res[j][bins[i][frame]] += float(l[j])
                    var[j][bins[i][frame]] += float(l[j])**2
                except KeyError :
                    res[j] = np.zeros(shape=nbins)
                    res[j][bins[i][frame]] += float(l[j])
                    var[j] = np.zeros(shape=nbins)
                    var[j][bins[i][frame]] += float(l[j])**2
        n += 1

# Set up out averages
avg = {}
vartotal = {}
for i in range(1,1+ncol) :
    avg[i] = 0.0
    vartotal[i] = 0.0

# Normalize the bin probability to the number of counts for each bin
# and multiply by the bin sums
totalprob = 0.0
for i in range(nbins) :
    totalprob += prob[i]
    if counts[i] > 0 :
        prob[i] /= counts[i]
invprob = 1.0 / totalprob

# Sum the bin totals * (normalized) bin probabilities
for i in range(nbins) :
    for j in range(1,1+ncol) :
        avg[j] += prob[i] * res[j][i]
        vartotal[j] += prob[i] * var[j][i] * invprob

# Print out the results
for i in range(1,1+ncol) :
    print "AVERAGE: %14.6f, STDEV: %14.6f, COLUMN: %5i"%(avg[i],np.sqrt(vartotal[i] - avg[i]**2),i)

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm
from math import log10
import sys
import os

class Read_H5 :
    def __init__ (self,filename,ndof) :
        try :
            self.file = h5.File(filename,"r")
            print "Loading %s"%filename
        except :
            print "Error opening %s"%filename
            sys.exit()
        self.ndof = ndof
        self.nexp = len(self.file['/Trajectories'].values())
        self.get_nconv()

    def get_dataset(self,dataset) :
        try :
            print "Reading %s"%dataset
            dat = self.file[dataset]
            results = []
            if len(dat.shape) > 1 :
                for i in range(len(dat[0])) :
                    results.append(dat[:,i])
            else :
                for i in range(len(dat)) :
                    results.append(dat[i])
            return results
        except KeyError :
            return -1

    def get_attribute(self,dataset,attribute) :
        dat = self.file[dataset]
        return dat.attrs.get(attribute)
    
    def get_nconv(self) :
        ensemble_dats = self.file['/Ensemble'].keys()
        nconv = 0
        for i in ensemble_dats :
            if "Conv-" in i :
                nconv += 1
        self.nconv = nconv
                
    def n_conv(self) :
        return self.nconv

    def conv_step_array(self) :
        step_xs = []
        steps = []
        step = 0
        for i in range(self.nconv) :
            set = '/Ensemble/Conv-%i'%i
            try :
                step = self.get_attribute('/Ensemble/Conv-%i'%i,'First frame read, last frame read')[1]
                steps.append(step*self.nexp)
            except KeyError :
                step = None
        step = self.get_attribute('/Ensemble','First frame read, last frame read')[1]
        steps.append(step*self.nexp)
        return np.array(steps)

    def conv_plot(self,ax,reference = None) :
        if reference == None :
            reference = self.get_dataset('/Ensemble/Probability')[-1]
        rsqr = []
        xs = self.conv_step_array()
        for i in range(self.nconv) :
            conv_i = np.array(self.get_dataset('/Ensemble/Conv-%i/Probability'%i)[-1])
            rsqr.append(np.sum((conv_i-reference)**2))
        rsqr.append(np.sum((reference-reference)**2))
        ax.plot(xs,rsqr)
        ax.set_yscale("log",nonposx='clip')
        ax.set_xlim(min(xs),max(xs)+1)
        ax.set_xscale("log",nonposx='clip')


def plot_2d(xs,ys,zs,ax,zmax=None,logscale=False,cb=True,zmin_avg=True) :
    if zmax == None :
        zmax = max(zs)
    minx,maxx = min(xs),max(xs)
    miny,maxy = min(ys),max(ys)
    xi = np.linspace(minx,maxx,360)
    yi = np.linspace(miny,maxy,360)
    if zmax == 100 or not logscale :
        zi = griddata(xs,ys,zs,xi,yi)
    else :
        Z = []
        zmin = 1e9
        for z in zs :
            try :
                Z.append(log10(z))
                if log10(z) < zmin :
                    zmin = int(np.floor(log10(z)))
            except ValueError :
                Z.append(-20)
            except :
                print "Unknown error trying to take the log of %s"%z
                sys.exit()
        zi = griddata(xs,ys,Z,xi,yi)
    myplot = ax.pcolorfast(xi,yi,zi,cmap='RdBu')
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    aspect_Ratio = (x1-x0)/(y1-y0)
    ax.set_aspect(aspect_Ratio)
    ax.set_xlim(-180,175)
    ax.set_ylim(-180,175)
    ax.set_xticks(range(-120,121,60))
    ax.set_yticks(range(-120,121,60))

    if zmax == 100 and not logscale :
        myplot.set_clim([0,100])
        if cb :
            cbar = plt.colorbar(myplot,shrink=0.70,format='%i',ticks=range(0,101,20))
            cbar.ax.set_yticklabels(['0%','20%','40%','60%','80%','>100%'])
    elif logscale :
        # Set the minimum value to the average magnitude
        if zmin_avg :
            zmin=int(np.floor(log10(np.average(zs))))
        myplot.set_clim(zmin,-1)
        if cb :
            cbar = plt.colorbar(myplot,format='%i',ticks=range(zmin,0,1))
            cblabels = []
            for i in range(zmin,0,1) :
                if i == zmin :
                    cblabels.append('<1E%i'%i)
                else :
                    cblabels.append('1E%i'%i)
            cbar.ax.set_yticklabels(cblabels)
    else :
        if cb :
            cbar = plt.colorbar(myplot,format='%1.E')
    return ax



import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm
from math import log10

"""
1) Read in hdf5 file
2) Extract probabilities
"""

def plot_2d(xs,ys,zs,ax,zmax = None) :
    if zmax == None : zmax = max(zs)
    minx,maxx=min(xs),max(xs)
    miny,maxy=min(ys),max(ys)
    xi = np.linspace(minx,maxx,360)
    yi = np.linspace(miny,maxy,360)
    zss = []
    if zmax == 100 :
        zi = griddata(xs,ys,zs,xi,yi)
    else :
        zmin = 1
        for z in zs:
            try :
                zss.append(log10(z))
                if z < zmin :
                    zmin = z
            except ValueError :
                zss.append(-12)
        zi = griddata(xs,ys,zss,xi,yi)
    myplot = ax.pcolorfast(xi,yi,zi,cmap='RdBu')
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect((x1-x0)/(y1-y0))
    ax.set_xlim(-180,175)
    ax.set_ylim(-180,175)
    ax.set_xticks(range(-120,121,60))
    ax.set_yticks(range(-120,121,60))

    if zmax == 100 :
        myplot.set_clim([0,zmax])
        cbar = plt.colorbar(myplot,shrink=0.70,format='%i',ticks=range(0,101,20))
        cbar.ax.set_yticklabels(['0%','20%','40%','60%','80%','>100%'])
    else :
        #cbar = plt.colorbar(myplot,shrink=0.70,format='%1.i')
        myplot.set_clim([-9,-1])
        cbar = plt.colorbar(myplot,shrink=0.70,format='%i',ticks=(-1,-3,-5,-7,-9))
        cbar.ax.set_yticklabels(['1E-1','1E-3','1E-5','1E-7','<1E-9'])


class read_hdf5 :
    def __init__ (self, filename, nwindows, ndof=1, doconv = True) :
        self.file = h5.File(filename,"r")
        self.ndof = ndof
        self.nwindows = nwindows
        self.nplots = 5
        self.prob1d = []
        self.prob2d = []
        self.rxs = []
        self.r = []
        self.nchunks = 0
        if doconv :
            if self.ndof == 1 :
                for i in range(len(self.file['/Ensemble'].values())) :
                    try :
                        self.prob1d.append(self.extract_prob_1D('Conv-%s/Probability'%i))
                        self.rxs.append( self.file['/Ensemble']['Conv-%i'%i].attrs.get('First frame read, last frame read')[1] * self.nwindows)
                        self.nchunks = i
                        print "Just read Conv-%i"%i
                    except:
                        break
            elif self.ndof == 2 :
                self.nwindows *= self.nwindows
                for i in range(len(self.file['/Ensemble'].values())) :
                    try :
                        self.prob2d.append(self.extract_prob_2D('Conv-%s/Probability'%i))
                        self.rxs.append( self.file['/Ensemble']['Conv-%i'%i].attrs.get('First frame read, last frame read')[1] * self.nwindows)
                        self.nchunks = i
                        print "Just read Conv-%i"%i
                    except :
                        break
            else :
                print "\nERROR!  Cannot handle %s degrees of freedom at this time"%self.ndof
                sys.exit()
        if self.ndof == 1 :
            print "Reading '/Ensemble/Probability' from %s"%filename
            self.prob1d.append(self.extract_prob_1D('Probability'))
            self.rxs.append( self.file['/Ensemble'].attrs.get('First frame read, last frame read')[1] * self.nwindows)
            self.nchunks += 1
        elif self.ndof == 2 :
            print "Reading '/Ensemble/Probability' from %s"%filename
            self.prob2d.append(self.extract_prob_2D('Probability'))
            self.rxs.append( self.file['/Ensemble'].attrs.get('First frame read, last frame read')[1] * self.nwindows)
            self.nchunks += 1
        else :
            print "\nERROR!  Cannot handle %s degrees of freedom at this time"%self.ndof
            sys.exit()

    def converge(self, analytic, ax) :
        if self.ndof == 1 : 
            dat = self.prob1d
        elif self.ndof == 2 : 
            dat = self.prob2d
        for i in range(self.nchunks+1) :
            try:
                self.r.append(self.cod(dat[i][self.ndof],analytic))
            except:
                break
        ax.plot(self.rxs,self.r)
        n = []
        for i in range(1,self.rxs[-1]) :
            n.append(1./(i))
        n2 = []
        for i in range(1,self.rxs[-1]) :
            n2.append(1./(i**.5))
        ax.plot(range(1,self.rxs[-1]),n,label='$1/n$')
        ax.plot(range(1,self.rxs[-1]),n2,label='$1/n^2$')
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        ax.set_xlabel('Number of frames over all windows')
        ax.set_ylabel(r'$\sum_i (\rho_{i_{WHAM}} - \rho_{i_{analytic}})^2$')
        ax.set_xlim(min(self.rxs),max(self.rxs)+1)
        ax.legend(loc='best')

    def plot_prob_1d(self,ax) : 
        dat = self.prob1d[-1]
        ax.plot(dat[0],dat[1])
        ax.set_xlim(-180,175)
        ax.set_xlabel('Angle (degrees)')
        ax.set_ylabel('Probability')
        ax.set_xticks(range(-120,121,60))
        
    def plot_prob_2d(self,ax,zmax = None) : 
        xs = self.prob2d[-1][0]
        ys = self.prob2d[-1][1]
        prob = np.array(self.prob2d[-1][2])
        plot_2d(xs,ys,prob,ax,zmax)

    def delta_prob_2d(self,ax,anaprob) : 
        xs = self.prob2d[-1][0]
        ys = self.prob2d[-1][1]
        
        dz = abs(anaprob - self.prob2d[-1][2])#/anaprob * 100.
        plot_2d(xs,ys,dz,ax)#,zmax=100)

    def cod(self, exp, ana):
        ssres = np.sum((exp-ana)**2)
        return ssres
        #sstot = np.sum((exp-np.average(exp))**2)
        #return 1 - ssres/sstot   
   
    def extract_prob_1D(self, path) : 
        dat = self.file['/Ensemble'][path]
        pdat = dat[()]
        prob = []
        x = []
        for i in range(len(pdat[:,[1]])) : 
            x.append(pdat[:,[0]][i][0])
            prob.append(pdat[:,[1]][i][0])
        x = np.array(x)
        prob = np.array(prob)
        return (x, prob)
        
    def extract_prob_2D(self, path) : 
        dat = self.file['/Ensemble'][path]
        pdat = dat[()]
        prob = []
        x = []
        y = []
        for i in range(len(pdat[:,[1]])) : 
            x.append(pdat[:,[0]][i][0])
            y.append(pdat[:,[1]][i][0])
            prob.append(pdat[:,[2]][i][0])
        x = np.array(x)
        y = np.array(y)
        prob = np.array(prob)
        return (x, y, prob)

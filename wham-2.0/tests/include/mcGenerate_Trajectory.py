#! /usr/bin/python

from numpy import array, arange, exp
from math import cos, sin, acos, atan2, pi
from pylab import random, randint, multinomial, shuffle, seed
import os
import sys
import matplotlib.pyplot as plt


class make_traj :
    def __init__(self, kb, T, fc, whamstep, dphi, windows, frames, potential, ws, we, datdir, stype = 'mc'):
        self.kb = kb
        self.T = T
        self.beta = 1/(kb*T)
        self.whamstep = whamstep
        self.dphi = dphi
        self.step = 1#whamstep
        self.windows = windows
        self.frames = frames
        self.potential = potential
        self.ws = ws
        self.we = we
        self.datdir = datdir
        self.R2D = 180/pi
        self.D2R = pi/180
        self.mk = 1
        self.fc = fc * self.D2R * self.D2R
        self.keep_every = 100
        if stype == 'mcmc' :
            self.mcmc = True
        else : 
            self.mcmc = False
        self.mc = not self.mcmc
        
    def mod2pi(self,phi):
        p = phi * self.D2R
        x = cos(p)
        y = sin(p)
        return atan2(y,x)*self.R2D
    
    def d_angle(self,phi1,phi2):
        phiprod = sin(phi1*self.D2R)*sin(phi2*self.D2R) + cos(phi1*self.D2R)*cos(phi2*self.D2R)
        if phiprod > 1 : 
            return 0
        return acos( phiprod )*self.R2D
    
    def vrb( self, phi ) :
        if self.potential == "dihedral" :
            # Proper dihedral potential
            V=1+cos(self.d_angle(3*phi,0))
            return V
        elif self.potential == "rb" :
            # Ryckaert-Bellemans alkane potential
            c=array((9.28,12.16,-13.12,-3.06,26.24,-31.5))
            Vrb=0
            for i in range(len(c)) :
                Vrb += c[i]*cos( self.d_angle( phi, 180. ) * self.D2R )**i
            return Vrb

    def bias_v(self, phi, phi0, fc, dphi=None):
        if dphi == None : dphi = self.dphi
        diffphi = self.d_angle(phi,phi0)
        if diffphi > dphi :
            ddp = self.d_angle(diffphi,dphi)
        else :
            ddp = 0
        return 0.5*fc*ddp*ddp
        
    def mod_v(self, phi, k=0):
        return k*2*sin(phi*self.D2R)   
    
    def bd_phi(self, phi, phi0, fc = None, step = None):
        if fc == None : fc = self.fc
        if step == None : step = self.step
        a = self.vrb(phi) + self.bias_v(phi,phi0,fc) + self.mod_v(phi,k=self.mk)
        b = self.vrb(phi+step) + self.bias_v(phi+step,phi0,fc) + self.mod_v(phi+step,k=self.mk)
        return 0.5 * step * ( exp(-self.beta*a) + exp(-self.beta*b) )
        
    def q_func(self, phi0, fc = None, step = None) :
        if fc == None : fc == self.fc
        if step == None : step = self.step
        totalV = 0
        for i in arange(-180,180,step) :
            totalV += self.bd_phi(i,phi0,fc,step)
        return totalV
    
    def generate_p_distr( self, x1phi0 = 0, x2phi0 = 0, fc = None ) :
        if fc == None : fc == self.fc
        xs1D = []
        prob = []
        pmf = []
        xs2D = []
        ys2D = []
        prob2D = []
        pmf2D = []
        q = self.q_func(x1phi0,fc,self.whamstep)
        for r in arange(-180,180,self.whamstep) :
            xs1D.append(r)
            prob.append( self.bd_phi(r,x1phi0,fc,self.whamstep)/q )
            pmf.append( self.vrb(r) + self.bias_v(r,x1phi0,fc) + self.mod_v(r,k=self.mk))
        for r1 in range(len(xs1D)) :
            x1_p = self.bd_phi(xs1D[r1],x1phi0,fc,self.self.whamstep)/q
            x1_pmf = self.vrb(xs1D[r1]) + self.bias_v(xs1D[r1],x1phi0,fc) + self.mod_v(xs1D[r1],k=self.mk)
            q2 = self.q_func(x2phi0,fc)
            for r2 in range(len(xs1D)) :
                xs2D.append(xs1D[r1])
                ys2D.append(xs1D[r2])
                x2_p = self.bd_phi(xs1D[r2],x2phi0,fc)/q2
                prob2D.append(x1_p * x2_p)
                x2_pmf = self.vrb(xs1D[r2]) + self.bias_v(xs1D[r2],x2phi0,fc) + self.mod_v(xs1D[r2],k=self.mk)
                pmf2D.append(x1_pmf + x2_pmf)
        return (array(xs1D), array(prob), array(pmf)),(array(xs2D),array(ys2D),array(prob2D),array(pmf2D))
 
    def unbiased_p_distr( self) :
        xs1D = []
        prob = []
        pmf = []
        xs2D = []
        ys2D = []
        prob2D = []
        pmf2D = []
        q = self.q_func(0,0,self.whamstep)
        for r in arange(-180,180,self.whamstep) :
            xs1D.append(r)
            prob.append( self.bd_phi(r,0,0,self.whamstep)/q )
            pmf.append( self.vrb(r) + self.bias_v(r,0,0) + self.mod_v(r,k=self.mk))
        for r1 in range(len(xs1D)) :
            for r2 in range(len(xs1D)) :
                xs2D.append(xs1D[r1])
                ys2D.append(xs1D[r2])
                prob2D.append(prob[r1] * prob[r2])
                pmf2D.append(pmf[r1] + pmf[r2])
        return (array(xs1D), array(prob), array(pmf)),(array(xs2D),array(ys2D),array(prob2D),array(pmf2D)) 
 
    def p_dist_1d(self, phi0 = 0, fc = None) :
        if fc == None : fc = self.fc
        q = self.q_func(phi0, fc)
        x = []
        y = []
        for i in arange(-180,180,self.step):
            x.append(i)
            y.append(self.bd_phi(i,phi0,fc)/q)
        #plt.plot(x,y,label=phi0)
        return array(y)
        
    def make_trajectory(self, distribution) : 
        if self.mcmc : 
            #print "Doing MCMC"
            return self.markov_trajectory(distribution)
        elif self.mc : 
            #print "Doing MC"
            return self.mc_trajectory(distribution)
        else : 
            print "Do not know what kind of sampling to do..."
            sys.exit()
            
    """ Markov chain Monte Carlo trajectory """
    def start_point(self, distribution) :
        x = multinomial(1,distribution)
        for i in range(len(x)) : 
            if x[i] == 1 :
                return i
                
    def markov_step(self,state,distribution):
        state = int(state)
        for i in range(self.keep_every) :
            if state == len(distribution) :
                state = 0
            if state == 0 :
                im = len(distribution)-1
            else :
                im = state-1
            if state == len(distribution) - 1:
                ip = 0
            else :
                ip = state+1   
            pi = distribution[state]
            pim = distribution[im]
            pip = distribution[ip]
            attempt = randint(0,2)
            if attempt == 0 : 
                if pim > pi : 
                    state = im
                else : 
                    z = random()
                    if z < pim/pi :
                        state = im
                    else :
                        state = state
            elif attempt == 1 :
                if pip > pi:
                    state = ip
                else : 
                    z = random()
                    if z < pip/pi :
                        state = ip
                    else :
                        state = state
        return state
        
    def markov_trajectory(self,distribution,state=None):
        if state == None :
            state = self.start_point(distribution)
        trj = []
        for i in range(self.frames):
            state = self.markov_step(state,distribution)
            angle = self.mod2pi((state+random())*self.step-180)
            trj.append(angle)
        return array(trj)
        
    """ Monte Carlo trajectory """
    def mc_trajectory(self, distribution) : 
        counts = multinomial(self.frames, distribution)
        samples = array( tuple( self.mc_counts2samples(counts) ) )
        shuffle(samples)
        return samples
        
    def mc_counts2samples( self, counts ) :
        values=[]
        bin=-180
        for value in counts:
            for i in arange(value) :
                values.append(bin+random()*self.step)
            bin +=  self.step
        return array(values)

    """ Generate outputs """
    def write_traj(self, phi0, traj, filelist, n, tag = "", dirnm = "") : 
        if not os.path.exists("%s/%s"%(self.datdir,dirnm)) : 
            os.makedirs("%s/%s"%(self.datdir,dirnm))
        filename = "%s/%s/%i%s.dat"%(self.datdir,dirnm,n,tag)
        datfile = open(filename, 'w')
        for i in traj : 
            datfile.write("%f\n"%i)
        datfile.close()
        filelist.write("%i %s %f %f %f %f\n"%(n,filename, phi0, self.dphi, self.fc * self.R2D * self.R2D, self.T))
        return filename
    
    def get_samples(self) : 
        if self.mcmc : 
            print "Doing MCMC"
        elif self.mc : 
            print "Doing MC"
        """ 1D Trajectory """
        filename1d = '1D_list.file'
        file1d = open(filename1d,'w')
        i =  0
        for w in arange(self.ws, self.we, self.windows) :
            print "Assigning probabilities %.2f"%w
            dist = self.p_dist_1d(phi0 = w, fc = self.fc)
            #plt.plot(arange(self.ws,self.we,self.step)-180,dist)
            traj = self.make_trajectory(dist)
            #plt.hist(traj,bins=72,normed=1,range=(-180,180),label="%i"%i)
            #plt.show()
            #plt.plot(traj,label=i)
            self.write_traj(w, traj, file1d, i,dirnm="dat_1d")
            i += 1
        """ 2D Trajectory """
        filename2d = '2D_list.file'
        file2d = open(filename2d, 'w')
        i = 0
        ii = 0
        for w1 in arange(self.ws, self.we, self.windows) : 
            j = 0
            dist1 = self.p_dist_1d(phi0 = w1, fc = self.fc)
            for w2 in arange(self.ws, self.we, self.windows) : 
                print "Assigning %.2f(%i), %.2f(%i)"%(w1,i,w2,j)
                dist2 = self.p_dist_1d(phi0 = w2, fc = self.fc)
                traj1 = self.make_trajectory(dist1)
                traj2 = self.make_trajectory(dist2)
                self.write_traj(w1, traj1, file2d, ii, tag = "-0d%i_%i"%(i,j), dirnm="dat_2d_%i"%i)
                self.write_traj(w2, traj2, file2d, ii, tag = "-1d%i_%i"%(i,j), dirnm="dat_2d_%i"%i)
                j += 1
                ii += 1
            i += 1
        return filename1d, filename2d
    
